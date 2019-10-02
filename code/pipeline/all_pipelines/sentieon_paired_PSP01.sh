#!/usr/bin/env bash

# *******************************************
# SENTIEON TUMOR/NORMAL VERSION 1.0 for paired analysis
# Dated 29-09-2019
# *******************************************
pipeline_version="PSP01"

# *******************************************
# ASSUMPTIONS
# two bam files are received from python submit scripts (somatic_setup.py)
# both bam files are produced by sentieon RNA / WGS pipeline and are already recalibrated
# *******************************************

# *******************************************
#                  SET-UP
# *******************************************

apps="/home/projects/HT2_leukngs/apps/github/code"

tumor_path=$1
normal_path=$2
destination=$3

# Number of threads
if [[ -z $4 ]]; then
        nt=28
        echo "# Running with default number of threads"
else
    echo Got threads $4
    nt=$4
fi

tumor_name=$(echo $(basename $tumor_path) | sed 's/.bam//g')
normal_name=$(echo $(basename $normal_path) | sed 's/.bam//g')
output_name=$(echo $(basename $destination))


echo -e "# Name \t Filepath: "
echo -e "# $tumor_name \t $tumor_path "
echo -e "# $normal_name \t $normal_path "
echo -e "# $output_name \t $destination "
echo "# Number of threads:" $nt
echo "# Reference version: b37"

set -x

mkdir $destination 

# Update with the location of the Sentieon software package
module load cbspythontools/1.0
SENTIEON_INSTALL_DIR=/services/tools/cbspythontools/1.0
export SENTIEON_LICENSE=localhost:8990

# Resources
reference_dir=/home/databases/gatk-legacy-bundles/b37
fasta=$reference_dir/human_g1k_v37_decoy.fasta
dbsnp=$reference_dir/dbsnp_138.b37.vcf
bedfile=/home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.filtered.bed

# *******************************************
#             PIPELINE
# *******************************************

## 1. Variant calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -r $fasta \
  -i $tumor_path \
  -i $normal_path \
  --algo TNscope \
  --tumor_sample $tumor_name --normal_sample $normal_name \
  --dbsnp $dbsnp \
  $destination/$output_name-TNscope.vcf.gz

echo -e "Finished the somatic variant calling" 

## 2. Statistics on all somatic variants
cd $destination 
$apps/ngs-tools/vcf_statistics.sh $output_name-TNscope.vcf.gz somatic

## 3. Analysis ...
vcftools --gzvcf $output_name-TNscope.vcf.gz --bed $bedfil --recode
nvariants=$(grep -v ^# -c out.recode.vcf)
echo "We get $nvariants somatic variants in the genes of interest"

# compress for VEP
echo "Compressing for VEP" 
bcftools view out.recode.vcf -Oz -o out.recode.vcf.gz
bcftools index out.recode.vcf.gz

# annotate
echo "Running VEP"
module purge    # pearl interference issues
$apps/ngs-tools/vep.sh out.recode
mv out.recode.vep.vcf.gz $output_name.vep.vcf.gz

# clean-up
echo "Deleting intermediate  files (out.recode*)" 
rm out.recode*
