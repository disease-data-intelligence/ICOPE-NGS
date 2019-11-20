#!/usr/bin/env bash

module load anaconda3/4.4.0
module load samtools/1.9
module load bedtools/2.28.0

echo "Using computerome modules samtools/1.9 and bedtools/2.28.0 and python3 modules listed below" 
bam=$1
sample=$(basename $bam | sed 's/.bam.*//')

start=`date +%s`

apps="/home/projects/HT2_leukngs/apps/github/code"
# make directory if not created already
destination="$sample".quality_reports
mkdir -p $destination

echo "# Getting coverage statistics for $sample in $destination"
bed=/home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.bed
genepanel=/home/projects/HT2_leukngs/data/references/general/300_genes_of_interest.txt


samtools bedcov $bed $bam >  $destination/cov.canonical.exons.bed
echo "# Summarizing exon and gene coverage"
$apps/quality/exon_coverage.py $destination/cov.canonical.exons.bed $genepanel
$apps/quality/gene_coverage.py $destination/cov.canonical.exons.bed $genepanel

intronbed=/home/projects/HT2_leukngs/data/references/hg37/UCSC_introns_hg37.bed
exonbed=/home/projects/HT2_leukngs/data/references/hg37/UCSC_exons_hg37.bed

samtools bedcov $intronbed $bam > $destination/introns.bed


echo "# Getting coverage pr. chromosome"
bedtools genomecov -ibam $bam -max 150 >  $destination/genome.cov
$apps/quality/chr_coverage.py  $destination/genome.cov 150

end=`date +%s`
runtime=$((end-start))

echo "# Runtime in seconds:" $runtime

