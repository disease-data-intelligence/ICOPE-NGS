#! /usr/bin/env bash
#!/bin/sh
# *******************************************
# Script to perform WGS variant calling
# using a single sample with fastq files
# This is VERSION 1.0 
# *******************************************
# Sentieon Germline version 1 
pipeline_version="PSG01"
# Full path to code repo 
apps="/home/projects/HT2_leukngs/apps/github/code"

# *******************************************
# ASSUMPTIONS
# fastq-files are name _R1 and R2 as the final part of the filename. After this only extension (can be zipped or not) 
# reads are aligned to hg37 decoy genome 
# *******************************************

set -x
#data_dir is  the directory of the input symlinks (symlinks should not be followed, done with -s option) 
data_dir="$( dirname "$(realpath -s $1)" )"
echo $data_dir
# use full path of input fastq 
fastq_1="$(realpath -s $1)"
fastq_2=$(sed 's/R1/R2/g' <<< "$fastq_1")


# Other settings
# add optional argument for number of threads, default is 28
if [[ -z $2 ]]; then
        nt=28
        echo "Threads: Running with defaults names"
else
    echo Got threads $2
    nt=$2
fi

echo "# Input files: $fastq_1 and $fastq_2"
echo "# Input number of threads:" $nt 

# Update with the location of the reference data files
fasta=/home/databases/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta						#ucsc.hg19_chr22.fasta
dbsnp=/home/databases/gatk-legacy-bundles/b37/dbsnp_138.b37.vcf							#dbsnp_135.hg19_chr22.vcf
known_1000G_indels=/home/databases/gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf
known_Mills_indels=/home/databases/gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf		#Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=localhost:8990

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/services/tools/cbspythontools/1.0

echo "SENTIEON INSTALL DIR="$SENTIEON_INSTALL_DIR 

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
# SENTIEON VERSION 1.0 code  
samplename=$(basename $fastq_1 | sed 's/.R1.*//')
sample="$samplename"."$pipeline_version"
group=$(zgrep -m 1 '@'  $fastq_1 | cut -d ':' -f3-4)
platform="ILLUMINA"


# ******************************************
# 0. Setup
# ******************************************
workdir=$data_dir/$sample
mkdir -p $workdir
logfile=$workdir/run.log
exec > $logfile 2>&1
echo "Copying script to" $workdir 
cp "$(readlink -f $0)" $workdir 
cd $workdir


# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $fasta -o sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam

# ******************************************
# 4. Indel realigner
# This step is optional for haplotyper-based caller like HC,
# but necessary for any pile-up based caller.
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels realigned.bam

# ******************************************
# 5. Base recalibration
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam  --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo ReadWriter recaled.bam

# ******************************************
# 5c. bam-statistics 
# ******************************************
# rename the final bam-file 
mv recaled.bam "$sample".bam 
mv recaled.bam.bai "$sample".bam.bai 

$apps/computerome/submit.py "$apps/ngs-tools/bam_statistics.sh "$sample".bam" --hours 15 -n "$sample"_bam_statistics -np 2 --no-numbering


# ******************************************
# 6. HC Variant caller
# Note: Starting GATK3.7, the default of emit_conf and call_conf is changed to 10.
# If you desire GATK-matching behavior, please change accordingly.
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=30 --call_conf=30 output-hc.vcf.gz


# ******************************************
# 7. Clean-up and VCF-stats submission  
# The final bam-file is the recaled.bam + *bai which may be used in other
# analysis 
# VCF-stats are found with vcf statistics and all quality is kept in
# quality_reports  
# ******************************************
# rename the final vcf-file to comply with naming scheme 
mv output-hc.vcf.gz "$sample"-hc.vcf.gz
mv output-hc.vcf.gz.tbi "$sample"-hc.vcf.gz.tbi
$apps/computerome/submit.py "$apps/ngs-tools/vcf_statistics.sh "$sample"-hc.vcf.gz" --name "$sample"_vcf_statistics -np 1 --no-numbering 

# remove all the files we don't want to keep: 
rm recal*
rm realigned.bam*
rm score.txt*
rm sorted.bam*
rm deduped.bam*

 
