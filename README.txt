#### INTRO
###########

## Supporting code for Heskett et al 2022: "Epigenetic Control of Hundreds of Chromosome-Associated lncRNA Genes Essential for Replication and Stability"
## All steps in the pipeline are based off of standard genomic processing and previously published software and analysis.

#### SOFTWARE TOOLS USED
########################
## Alignment: BWA 0.7.17-r1188, STAR 020201
## Processing alignments: SAMTools 1.9, BCFtools 1.9, BEDtools v2.27.1
## Statistics and plotting: Python 3.6.7, Pandas 0.24.2, NumPy 1.16.2, MatPlotLib 3.1.0

#### GENERAL DESCRIPTION OF PIPELINE
####################################
## Allele-specific RNA-seq analysis of ASAR lncRNAs and coding genes is achieved by comparing strand-specific
## nuclear-enriched, totalRNA with Ribo-minus, RNA-seq libraries to a fully haplotype-resolved reference genotype
## to enumerate allele-specific reads. Alignment is performed with STAR, allele-specific read counting with
## BCFtools mpileup, and the binomial test and plotting is done with Python 3. 
## 
## 1. Align RNA-seq reads for each sample using STAR, remove duplicates, and filter low quality reads with SAMtools.
## 2. Identify strand-specific Very Long Non-Coding RNAs (VLINCs) or non-coding "Transcribed Loci" for each sample or group of samples. 
##    Identifying vlincs is originally described in St. Laurent et al 2013, Genome Biology. Also see "Methods: VlincRNA identification
##    and  expression profiling" in M. Caron et al 2018, PLoS one.
## 3. Obtain strand specific bed file with coding genes of interest.
## 4. Run BCFtools mpileup using the haplotype-resolved VCF file to select all sites with heterozygous SNPs. 
## 5. Use BEDtools to overlap both the vlincRNA file and the coding genes file with the BCFtools mpileup output
##    to get strand specific read counts for all vlincRNAs and coding genes. 
## 6. Use Python 3 to get FDR corrected binomial-test q-values and make plots. 
##
## Allele-specific Repli-seq analysis is achieved by comparing repli-seq reads to a fully haplotype-resolved
## reference genotype to enumerate allele-specific reads in sliding genomic windows across the genome. Replication timing
## analysis is desribed in Marchal et al. 2019, Nature protocols.
##
## 1. Align early and late fraction repli-seq reads with BWA, remove duplicates, and filter low quality reads with SAMtools.
## 2. Count allele-specific reads using BCFtools mpileup using the haplotype-resolved VCF file to select all sites
##    with heterozygous SNPs.  
## 3. Generate sliding genomic windows with BEDtools makewindows and the reference genome fasta.fai index file. 
## 4. Use BEDtools with the count option to count allele-specific reads in each sliding window. RT profiles are
##    defined as the log2 Ratio of Early to the Late fraction for each haplotype. 
## 5. Use Python 3 to calculate the sandard deviation of RT within each window and make plots. 

#### SOFTWARE METHODS PAPERS CITED
## Identifying VLINCs or "TLs": St. Laurent et al 2013, Genome Biology, M. Caron et al 2018, PLoS one.
## Repli-Seq Analysis: Marchal et al. 2019, Nature protocols.
## Standard Genomics processing: Bedtools: Quinlan et al Bioinformatics 2010. SAMtools: Heng Li et al Bioinformatics, 2009. 
## BWA: Heng Li et al Bioinformatics 2009. Allele-specific analysis: S. Castel et al 2015 Genome Biology. 
