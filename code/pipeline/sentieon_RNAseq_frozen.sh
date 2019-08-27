#!/bin/sh
# *******************************************
# Script to perform RNA seq variant calling
# *******************************************


# SENTIEON TUMOR VERSION 1.0 for RNA seq
pipeline_version="PST01"

# Full path to code repo
apps="/home/projects/HT2_leukngs/apps/github/code"


# *******************************************
# ASSUMPTIONS
# fastq-files are name _R1 and R2 as the final part of the filename. After this only extension (can be zipped or not)
# reads are aligned to hg37 decoy genome
# *******************************************

set -x
#data_dir is  the directory of the input symlinks (symlinks should not be followed,which is configured with -s option)
data_dir="$( dirname "$(realpath -s $1)" )"
fastq_1="$(realpath -s $1)" 
fastq_2=$(sed 's/R1/R2/g' <<< "$fastq_1")

# other settings
nt=$2 #number of threads to use in computation

echo "# Input files: $fastq_1 and $fastq_2"
echo "# Input number of threads:" $nt
echo "# Aligning to b37"

reference_dir=/home/databases/gatk-legacy-bundles/b37

# Update with the location of the reference data files
fasta=$reference_dir/human_g1k_v37_decoy.fasta								#ucsc.hg19_chr22.fasta
dbsnp=$reference_dir/dbsnp_138.b37.vcf										#dbsnp_135.hg19_chr22.vcf
known_1000G_indels=$reference_dir/1000G_phase1.indels.b37.vcf			#1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf
known_Mills_indels=$reference_dir/Mills_and_1000G_gold_standard.indels.b37.vcf			#Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf

# Update with the location of the Sentieon software package
module load cbspythontools/1.0
SENTIEON_INSTALL_DIR=/services/tools/cbspythontools/1.0
export SENTIEON_LICENSE=localhost:8990

# RNA specific resources
module load star/2.7.0c
echo "Loaded star version 2.7.0c" 
genomeDir=/home/projects/HT2_leukngs/data/references/hg37/star_genome_hg37

# It is important to assign meaningful names in actual cases meaning that we will follow the usual naming scheme! 
# It is particularly important to assign different read group names.

samplename=$(basename $fastq_1 | sed 's/.R1.*//')
sample="$samplename"."$pipeline_version"
RG=$(zgrep -m 1 '@'  $fastq_1 | cut -d ':' -f3-4)
platform="ILLUMINA"
platformUnit=$(zgrep -m 1 '@'  $fastq_1 | cut -d ':' -f3)

# 0. Setup
# ******************************************
workdir=$data_dir/$sample
mkdir -p $workdir
logfile=$workdir/run.log
exec > $logfile 2>&1
echo "Copying executed script to" $workdir
cp "$(readlink -f $0)" $workdir
cd $workdir

# ******************************************
#$ 1. Mapping with STAR and sorting
# Sorting happens by coordinate
# We use 'basic' two-pass mode which means that STAR will perform the 1st pass mapping, then it will automatically extract junctions, insert them into the genome index, and,
# finally, re-map all reads in the 2nd mapping pass.
# outSAMmapqUnique makes sure that the mapping qualities (255) for uniquely mapped reads is reset to 60
# ******************************************

STAR --genomeDir $genomeDir \
  --readFilesIn $fastq_1 $fastq_2 \
  --readFilesCommand zcat \
  --runThreadN $nt \
  --twopassMode Basic \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix sorted"." \
  --outSAMmapqUnique 60

# ******************************************
# 1b. Sorting
# outSAMtype does not work, so we use samtools for sorting by coordinate
# ******************************************
module load samtools/1.9
echo "Loaded samtools 1.9 for sorting, indexing and adding read group" 
samtools addreplacerg sorted.Aligned.out.bam -r ID:$RG -r PU:$platformUnit -r SM:$sample -r PL:$platform -m overwrite_all -o rg.bam
samtools sort rg.bam -o sorted.bam --threads $nt
samtools index sorted.bam


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
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam


# ******************************************
# 4. Split reads at junction
# This step splits the RNA reads into exon segments by getting rid of Ns while maintaining grouping information, and hard-clips any sequences overhanging into the intron regions.
# Additionally, the step will reassign the mapping qualities from STAR to be consistent with what is expected in subsequent steps by converting from quality 255 to 60.
# ^ this is not really necessary (this is from the website / manual)
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --algo RNASplitReadsAtJunction --reassign_mapq 255:60 splitted.bam

# ******************************************
# 5. Indel realigner
# This step is optional for haplotyper-based caller like HC,
# but necessary for any pile-up based caller.
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i splitted.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels realigned.bam

# ******************************************
# 6. Base recalibration
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 6b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo ReadWriter recaled.bam

# ******************************************
# 6c. bam-statistics
# ******************************************
# rename the final bam-file
mv recaled.bam "$sample".bam
mv recaled.bam.bai "$sample".bam.bai

$apps/computerome/submit.py "$apps/ngs-tools/bam_statistics.sh "$sample".bam" --hours 15 -n "$sample"_bam_statistics -np 2 --no-numbering


# ******************************************
# 6a. HC Variant caller
# Note: Starting GATK3.7, the default of emit_conf and call_conf is changed to 10.
# If you desire GATK-matching behavior, please change accordingly.
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper --trim_soft_clip --call_conf 20 --emit_conf 20 -d $dbsnp output-hc.vcf.gz


# ******************************************
# 6b. Variant calling with TNscope without normal sample
# Here we could include a 'panel of normals' --pon (VCF-file) which should substitute the normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo TNscope --tumor_sample $sample -d $dbsnp --trim_soft_clip --disable_detector sv output-TNScope.vcf.gz


# ******************************************
# 7. Clean-up and VCF-stats submission
# The final bam-file is the recaled.bam + *bai which may be used in other
# analysis
# VCF-stats are found with vcf statistics and all quality is kept in
# quality_reports
# ******************************************
# rename the final vcf-file to comply with naming scheme
mv output-hc.vcf.gz "$sample"-hc.vcf.gz
mv output-TNScope.vcf.gz "$sample"-TNScope.vcf.gz
mv output-hc.vcf.gz.tbi "$sample"-hc.vcf.gz.tbi
mv output-TNScope.vcf.gz.tbi "$sample"-TNScope.vcf.gz.tbi 


$apps/computerome/submit.py "$apps/ngs-tools/vcf_statistics.sh "$sample"-hc.vcf.gz hc" --name "$sample"-hc-vcf_statistics -np 1 --no-numbering

# remove all the files we don't want to keep:
rm recal* 
rm splitted.bam*
rm realigned.bam*
rm score.txt*
rm sorted.Aligned*
rm rg.bam*
rm deduped.bam*
rm sorted.bam* 

