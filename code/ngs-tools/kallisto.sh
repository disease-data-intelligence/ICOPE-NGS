#! /usr/bin/env bash

module load kallisto/0.46.0
index=/home/projects/HT2_leukngs/data/references/homo_sapiens/transcriptome.idx
fastq_1=$1
fastq_2=$(sed 's/R1/R2/g' <<< "$fastq_1")

kallisto quant -i $index -o kallisto --fusion $fastq_1 $fastq_2 
