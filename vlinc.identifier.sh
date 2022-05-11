#!/bin/bash

input=$1 ## bam file
out_dir=$2 ## output directory
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension

library_size=$(samtools view -c -F 260 $input) # get library size
echo $library_size

## set merging parameters
first_merge=$3
second_merge=$4
filter_size=$5

## get plus strand
samtools view -b -f 131 -F 16 $input > $out_dir$filename.fwd1.bam
samtools index $out_dir$filename.fwd1.bam

# now get first in pair mapping to reverse strand
samtools view -b -f 83 $input > $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.fwd2.bam

# now combine, this should contail all plus strand genes
samtools merge -f $out_dir$filename.plus.bam $out_dir$filename.fwd1.bam $out_dir$filename.fwd2.bam
samtools index $out_dir$filename.plus.bam

## get minus strand as same as above
samtools view -b -f 147 $input > $out_dir$filename.rev1.bam
samtools index $out_dir$filename.rev1.bam

samtools view -b -f 67 -F 16 $input > $out_dir$filename.rev2.bam
samtools index $out_dir$filename.rev2.bam

samtools merge -f $out_dir$filename.minus.bam $out_dir$filename.rev1.bam $out_dir$filename.rev2.bam
samtools index $out_dir$filename.minus.bam

## now do coverage separately for each strand
bedtools genomecov -bg -split -ibam $out_dir$filename.plus.bam > $out_dir$filename.plus.cov.bed
bedtools genomecov -bg -split -ibam $out_dir$filename.minus.bam > $out_dir$filename.minus.cov.bed

## Remove all overlap with known exons and introns for coding genes (whole gene) in a strand specific manner. Remove poor mapping regions

bedtools subtract -a $out_dir$filename.plus.cov.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/encode.blacklisted.nochr.hg19.bed \
  | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.plus.nochr.hg19.bed \
  > $out_dir$filename.plus.subtract.whole.gene.bed

bedtools subtract -a $out_dir$filename.minus.cov.bed -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/encode.blacklisted.nochr.hg19.bed \
  | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.minus.nochr.hg19.bed \
  > $out_dir$filename.minus.subtract.whole.gene.bed

## merge segments that are separated by --first_merge.Then subtract the cds again incase you merged over some tiny coding exons.
bedtools merge -i $out_dir$filename.plus.subtract.whole.gene.bed -d $first_merge | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.plus.nochr.hg19.bed \
  > $out_dir$filename.plus.subtract.merge.$first_merge.bed
bedtools merge -i $out_dir$filename.minus.subtract.whole.gene.bed -d $first_merge | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.minus.nochr.hg19.bed \
  > $out_dir$filename.minus.subtract.merge.$first_merge.bed

# now merge segments if theyre less than $second_merge
bedtools merge -i $out_dir$filename.plus.subtract.merge.$first_merge.bed -d $second_merge | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.plus.nochr.hg19.bed \
  > $out_dir$filename.plus.$first_merge.$second_merge.bed
bedtools merge -i $out_dir$filename.minus.subtract.merge.$first_merge.bed -d $second_merge | bedtools subtract -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.minus.nochr.hg19.bed \
   > $out_dir$filename.minus.$first_merge.$second_merge.bed

## bedtools coverage with the strand option (-s) doesnt work as I would expect it to, so I am avoiding this. 
## Filters by min size 50kb, then removes anything that has -f amount of overlap (0.25) with any codig exon. 
## Get read coverage again from the BAM
## Overlaps LINE elements that have 90% coverage of the LINE
## Basic calculation of RPKM and print to file
## this step requires a large amount of memory
awk -v var="$filter_size" -v var2="$filename" '{OFS="\t"} $3-$2>=var{print $1, $2, $3, i++"_plus_"var2}' $out_dir$filename.plus.$first_merge.$second_merge.bed |
	bedtools intersect -f 0.25 -v -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.nochr.hg19.bed |
	bedtools coverage -c -F 0.5 -a stdin -b $out_dir$filename.plus.bam |
	bedtools coverage -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed -F 0.9 |
	awk -v var3="$library_size" 'OFS="\t"{print $1, $2, $3, $4, $5 / ( ($3-$2)/1000 ), "+", $9}' > $out_dir$filename.intergenic.plus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed

awk -v var="$filter_size" -v var2="$filename" -v var3="$library_size" '{OFS="\t"} $3-$2>=var{print $1, $2, $3, i++"_minus_"var2}' $out_dir$filename.minus.$first_merge.$second_merge.bed |
	bedtools intersect -f 0.25 -v -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.known.genes.cds.nochr.hg19.bed |
	bedtools coverage -c -F 0.5 -a stdin -b $out_dir$filename.minus.bam |
	bedtools coverage -a stdin -b /home/groups/Spellmandata/heskett/replication.rnaseq/annotation.files/ucsc.L1.filtered.hg19.bed -F 0.9 |
	awk -v var3="$library_size" 'OFS="\t"{print $1, $2, $3, $4, $5 / ( ($3-$2)/1000 ), "-", $9}' > $out_dir$filename.intergenic.minus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed

## create a	new file specifically in the genome browser format
awk 'OFS="\t"{print "chr"$1, $2, $3, $4, 0, $6}' $out_dir$filename.intergenic.plus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI > $out_dir$filename.intergenic.$first_merge.$second_merge.$filter_size.vlinc.discovery.plus.browser.bed
awk 'OFS="\t"{print "chr"$1, $2, $3, $4, 0, $6}' $out_dir$filename.intergenic.minus.$first_merge.$second_merge.$filter_size.vlinc.discovery.bed | grep -v chrGL | grep -v chrKI > $out_dir$filename.intergenic.$first_merge.$second_merge.$filter_size.vlinc.discovery.minus.browser.bed


### remove intermediate files to avoid taking up tons of space. Intermediate files are required to reduce
### the amount of memory needed. 
rm $out_dir$filename.plus.cov.bed
rm $out_dir$filename.minus.cov.bed
rm $out_dir$filename.plus.subtract.whole.gene.bed
rm $out_dir$filename.minus.subtract.whole.gene.bed
rm $out_dir$filename.plus.subtract.merge.$first_merge.bed
rm $out_dir$filename.minus.subtract.merge.$first_merge.bed
rm $out_dir$filename.plus.$first_merge.$second_merge.bed
rm $out_dir$filename.minus.$first_merge.$second_merge.bed
rm $out_dir$filename.fwd1.bam*
rm $out_dir$filename.fwd2.bam*
rm $out_dir$filename.rev1.bam*
rm $out_dir$filename.rev2.bam*
