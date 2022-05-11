#!/bin/bash

fq1=$1
fq2=$2
name=$3
ref='/fasta/ref'
out_dir=$4
readgroup='@RG\tID:'$name'\tSM:'1'\tPL:'illumina'\tLB:'$name
bwa mem -R $readgroup -t 8 -Ma $ref $fq1 $fq2 > $out_dir$name.sam
