## This script takes a genes or ASAR bed file, and intersects with allele-specific read counts at heterozgyous SNP sites
## to get the allele-specific read counts for each gene or ASAR

bedtools map -a ucsc.genes.plus.bed -b sample.allele.counts.plus.bed \
   -o sum,sum -c 6,7 > sample.protein.coding.plus.rna.bed
bedtools map -a ucsc.genes.minus.bed -b sample.allele.counts.minus.bed \
   -o sum,sum -c 6,7 > sample.protein.coding.minus.rna.bed
cat sample.allele.counts.plus.bed sample.allele.counts.minus.bed | \
     sort -k1,1 -k2,2n | grep -Fv "." > sample.protein.coding.all.counts.bed
