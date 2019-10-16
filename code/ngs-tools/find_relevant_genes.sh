#!/usr/bin/env bash

# *****************************************************************
# Script to go through variant filtering on interesting genes
# *****************************************************************

# ASSUMPTIONS: Input is .vcf og .vcf.gz

# bedfile with interesting genes in ALL
bedfile=/home/projects/HT2_leukngs/data/references/hg37/collected_hg37_ALL_genes.bed
# bedfile with pathways genes:
# bedfile_pathways

input=$1
sample=$(basename $input | sed 's/.vcf.*//')

echo -e "Input \t Sample-name"
echo -e "$input \t $sample"

# subset list of variants with our gene list (bed-file)
module load vcftools/0.1.16
echo "Loaded vcftools/0.1.16"
vcftools --gzvcf $sample".vcf.gz" --bed $bedfile --recode
nvariants=$(grep -v ^# -c out.recode.vcf)
echo "We get $nvariants variants in the genes of interest (all ALL genes)"

# compress for VEP
echo "Compressing for VEP"
module load bcftools/1.9
echo "Loaded bcftools/1.9"
bcftools view out.recode.vcf -Oz -o out.recode.vcf.gz
bcftools index out.recode.vcf.gz

# annotate
echo "Running VEP"
module purge    # pearl interference issues
# some people use apps alias in other settings so it is overwritten again
apps="/home/projects/HT2_leukngs/apps/github/code"
$apps/ngs-tools/vep.sh out.recode.vcf.gz
$apps/ngs-tools/filter_vep.sh out.recode.vep.vcf.gz

module load tools
module load bcftools/1.9
echo "Loaded bcftools/1.9"
echo "Output stored as $sample.filter.vep.vcf.gz"
bcftools view out.recode.vep.filter.vcf -Oz -o out.recode.vep.filter.vcf.gz
mv out.recode.vep.filter.vcf.gz  $sample.filter.vep.vcf.gz
echo "Creating index ..." 
bcftools index $sample.filter.vep.vcf.gz 

module load anaconda3/4.4.0
$apps/ngs-tools/summarize_vep_variants.py -samples $sample.filter.vep.vcf.gz -outfile $sample


# Cleaning up 
set -x 
rm out*
echo "## FINISHED FIRST PART, NOW LOOKING ONLY AT PATHWAY GENES ## "

bedfile=/home/projects/HT2_leukngs/data/references/hg37/hg37_pathways.bed
# subset list of variants with our gene list (bed-file)
module load vcftools/0.1.16
echo "Loaded vcftools/0.1.16"
vcftools --gzvcf $sample".vcf.gz" --bed $bedfile --recode
nvariants=$(grep -v ^# -c out.recode.vcf)
echo "We get $nvariants variants in the genes of interest (all pathways ALL genes)"

# compress for VEP
echo "Compressing for VEP"
module load bcftools/1.9
echo "Loaded bcftools/1.9"
bcftools view out.recode.vcf -Oz -o out.recode.vcf.gz
bcftools index out.recode.vcf.gz

# annotate
echo "Running VEP"
module purge    # pearl interference issues
# some people use apps alias in other settings so it is overwritten again
apps="/home/projects/HT2_leukngs/apps/github/code"
$apps/ngs-tools/vep.sh out.recode.vcf.gz
$apps/ngs-tools/filter_vep.sh out.recode.vep.vcf.gz

module load tools
module load bcftools/1.9
echo "Loaded bcftools/1.9"
echo "Output stored as $sample.pathways.filter.vep.vcf.gz"
bcftools view out.recode.vep.filter.vcf -Oz -o out.recode.vep.filter.vcf.gz
mv out.recode.vep.filter.vcf.gz  $sample.pathways.filter.vep.vcf.gz
echo "Creating index ..."
bcftools index $sample.pathways.filter.vep.vcf.gz

module load anaconda3/4.4.0
$apps/ngs-tools/summarize_vep_variants.py -samples $sample.pathways.filter.vep.vcf.gz -outfile $sample