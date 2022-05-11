#!/bin/bash

## start from sorted, index and rmduped BAM file.

input=$1
b=$(basename $input) ## removes /path/to/file
filename=${b%.*} ### removes file extension
out_dir=$2

# Phased genotype bed file needed: PHASED.GENOTYPE.HET.SITES.bed
# genome reference file needed
# allele specific read counts are output into a VCF

bcftools mpileup -R /home/groups/Spellmandata/heskett/replication.rnaseq/scripts/PHASED.GENOTYPE.HET.SITES.bed \
  -f reference.file.fasta \
  -a DP,AD \
  -q 20 \
  $1 > $out_dir$filename.allele.counts.vcf

# get allele specific counts of all het sites into a table
gatk VariantsToTable -V $out_dir$filename.allele.counts.vcf -O $out_dir$filename.table -F CHROM -F POS -F REF -F ALT -GF AD

## turn this table into a bed file so that it can be intersected with the haplotype resolved file and used with BEDtools
tail -n +2 $out_dir$filename.table | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' | 
bedtools intersect -a stdin -b PHASED.GENOTYPE.HET.SITES.bed -wa -wb > $out_dir$filename.allele.counts.bed
