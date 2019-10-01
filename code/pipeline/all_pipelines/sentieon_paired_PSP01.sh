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
echo -e "# $tumor_name" \t "$tumor_path "
echo -e "# $normal_name" \t "$normal_path "
echo -e "# $output_name" \t "$destination " 
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
  --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
  $destination/$output_name-TNscope.vcf.gz

## 2. Statistics
cd $destination 
$apps/ngs-tools/vcf_statistics.sh $output_name-TNscope.vcf.gz somatic

## 3. Analysis ...