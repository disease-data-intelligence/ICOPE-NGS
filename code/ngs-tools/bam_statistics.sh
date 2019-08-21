#!/usr/bin/env bash

module load anaconda3/4.4.0
module load samtools/1.9
module load bedtools/2.28.0

echo "Using computerome modules samtools/1.9 and bedtools/2.28.0 and python3 modules listed below" 
bam=$1

start=`date +%s`

apps="/home/projects/HT2_leukngs/apps/github/code"
mkdir -p quality_reports

echo "# Getting coverage statistics"
bed=/home/projects/HT2_leukngs/data/references/ hg37/USCS.hg37.canonical.exons.filtered.bed
genepanel=/home/projects/HT2_leukngs/data/references/general/genepanel.txt

samtools bedcov $bed $bam > quality_reports/cov.canonical.exons.bed
echo "# Summarizing exon and gene coverage"
$apps/quality/exon_coverage.py quality_reports/cov.canonical.exons.bed $genepanel 
$apps/quality/gene_coverage.py quality_reports/cov.canonical.exons.bed $genepanel

echo "# Getting coverage pr. chromosome"
bedtools genomecov -ibam $bam -max 150 > quality_reports/genome.cov
$apps/quality/chr_coverage.py quality_reports/genome.cov

end=`date +%s`
runtime=$((end-start))

echo "# Runtime in seconds:" $runtime
