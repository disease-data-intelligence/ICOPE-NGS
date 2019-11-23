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
intronbed=/home/projects/HT2_leukngs/data/references/hg37/UCSC_introns_hg37.bed
genepanel=/home/projects/HT2_leukngs/data/references/general/315_genes_of_interest.txt
repair_genes=/home/projects/HT2_leukngs/data/references/general/DNA_repair_genes_core.txt

echo "# Summarizing exon and gene coverage"
start_bedcov=`date +%s`
samtools bedcov $bed $bam > $destination/cov.canonical.exons.bed
end_bedcov=`date +%s`
runtime=$((end_bedcov-start_bedcov))
echo "Finished bedov in $runtime"

$apps/quality/combined_coverage.py -in $destination/cov.canonical.exons.bed -panel $genepanel -out "genes_of_interest"
$apps/quality/combined_coverage.py -in $destination/cov.canonical.exons.bed -panel $repair_genes -out "repair_genes"

# $apps/quality/exon_coverage.py $destination/cov.canonical.exons.bed $genepanel
# $apps/quality/gene_coverage.py $destination/cov.canonical.exons.bed $genepanel

# samtools bedcov $intronbed $bam > $destination/cov.introns.bed
$apps/quality/combined_coverage.py -in $bam -bed $intronbed -intron

echo "# Getting coverage pr. chromosome"
#bedtools genomecov -ibam $bam -max 150 >  $destination/genome.cov
$apps/quality/chr_coverage.py  $bam 150

end=`date +%s`
runtime=$((end-start))

echo "# Runtime in seconds:" $runtime

