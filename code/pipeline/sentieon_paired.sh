#!/usr/bin/env bash




# *******************************************
# ASSUMPTIONS
# two bam files are received from python submit scripts
# both bam files are produced by sentieon RNA / WGS pipeline and already recalibrated
# *******************************************

# *******************************************
#                  SET-UP
# *******************************************

# SENTIEON TUMOR/NORMAL VERSION 1.0 for paired analysis
# Dated 29-09-2019
pipeline_version="PSP01"

apps="/home/projects/HT2_leukngs/apps/github/code"

tumor_path=$1
germline_path=$2
destination=$3

# Number of threads
if [[ -z $4 ]]; then
        nt=5
        echo "# Running with default number of threads"
else
    echo Got threads $4
    nt=$4
fi

tumor_name=$(echo $(basename $tumor_path) | sed 's/.bam//g')
normal_name=$(echo $(basename $normal_path) | sed 's/.bam//g')
output_name="$tumor_name"$normal_name"$pipeline_version"


echo "# Input files: $tumor_path and $germline_path"
echo "# Names: $tumor_name, $normal_name, $output_name"
echo "# Results are stored in $destination"
echo "# Number of threads:" $nt
echo "# Reference version: b37"

set -x
#data_dir is  the directory of the input symlinks (symlinks should not be followed,which is configured with -s option)
data_dir="$( dirname "$(realpath -s $1)" )"

# Update with the location of the Sentieon software package
module load cbspythontools/1.0
SENTIEON_INSTALL_DIR=/services/tools/cbspythontools/1.0
export SENTIEON_LICENSE=localhost:8990

# Resources
reference_dir=/home/databases/gatk-legacy-bundles/b37
fasta=$reference_dir/human_g1k_v37_decoy.fasta								#ucsc.hg19_chr22.fasta
dbsnp=$reference_dir/dbsnp_138.b37.vcf

# *******************************************
#             PIPELINE
# *******************************************

## 1. Variant calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -r $fasta \
  -i $tumor_path \
  -i $germline_path \
  --algo TNscope \
  --tumor_sample $tumor_name --normal_sample $normal_name \
  --dbsnp $dbsnp \
  $destination/$output_name-TNscope.vcf.gz

## 2. Statistics
$apps/ngs-tools/vcf_statistics.sh $destination/$output_name-TNscope.vcf.gz somatic

## 3. Analysis ...
