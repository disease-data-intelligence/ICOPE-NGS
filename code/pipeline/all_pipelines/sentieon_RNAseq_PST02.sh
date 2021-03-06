#!/bin/sh
# *******************************************
# Script to perform RNA seq variant calling
# *******************************************


# SENTIEON TUMOR VERSION 2.0 for RNA seq
# Dated 29-09-2019
pipeline_version="PST02"

# Full path to code repo
apps="/home/projects/HT2_leukngs/apps/github/code"


# *******************************************
# ASSUMPTIONS
# fastq-files are named R1 and R2 as the final part of the filename. After this only extension (can be zipped or not)
# reads are aligned to hg37 decoy genome
# *******************************************

set -x
#data_dir is  the directory of the input symlinks (symlinks should not be followed, which is configured with -s option)
data_dir="$( dirname "$(realpath -s $1)" )"
fastq_1="$(realpath -s $1)" 
fastq_2=$(sed 's/R1/R2/g' <<< "$fastq_1")

# other settings
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
echo "# Aligning to b37"

reference_dir=/home/databases/gatk-legacy-bundles/b37

# Update with the location of the reference data files
fasta=$reference_dir/human_g1k_v37_decoy.fasta
dbsnp=$reference_dir/dbsnp_138.b37.vcf
known_1000G_indels=$reference_dir/1000G_phase1.indels.b37.vcf
known_Mills_indels=$reference_dir/Mills_and_1000G_gold_standard.indels.b37.vcf

# Update with the location of the Sentieon software package
module load cbspythontools/1.0
SENTIEON_INSTALL_DIR=/services/tools/cbspythontools/1.0
export SENTIEON_LICENSE=localhost:8990

# RNA specific resources
module load gcc/8.2.0
module load star/2.7.2b
echo "Loaded gcc/8.2.0 and star version 2.7.2b"
gtf=/home/projects/HT2_leukngs/data/references/hg37/Homo_sapiens.GRCh37.75_filtered.gtf
genomeDir=/home/projects/HT2_leukngs/data/references/hg37/star_genome_2.7.2b

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
script_name=$(basename $(readlink -f $0))
cp "$(readlink -f $0)" $workdir/"$sample".$script_name
cd $workdir

# ******************************************
# 0a. Checking input data
# ******************************************
nr_lines1=$(zcat $fastq_1 | wc -l)
nr_lines2=$(zcat $fastq_2 | wc -l)
raw_reads1=$(echo "scale=2 ; $nr_lines1 / 4" | bc)
raw_reads2=$(echo "scale=2 ; $nr_lines2 / 4" | bc)
raw_read_threshold_soft=20000000
raw_read_threshold_hard=10000000

echo $raw_reads1 and $raw_reads2  ...
if [ $raw_reads1 -lt $raw_read_threshold_soft ]; then
    echo "Warning! The number of reads might not be sufficient ... !"
    if [ $nr_lines1 -lt $raw_read_threshold_hard ]; then
    echo "Number of reads is insufficient, exiting ..."
    exit
    fi
fi
if [ $raw_reads1 == $raw_reads2 ]; then
    echo "Number of reads for each sample are equal"
    else
    echo "Warning: Number of reads are not equivalent, fastqc files might be corrupted"
fi

echo "Continuing with alignment ... "


# ******************************************
# 1a. Mapping with STAR and sorting
# Sorting happens by coordinate
# We use 'basic' two-pass mode which means that STAR will perform the 1st pass mapping, then it will automatically
# extract junctions, insert them into the genome index, and, finally, re-map all reads in the 2nd mapping pass.
# outSAMmapqUnique makes sure that the mapping qualities (255) for uniquely mapped reads is reset to 60 to match bwa output
# we name output files according to the sample name as the directories should have sample naming as well
# ******************************************

STAR --genomeDir $genomeDir \
  --readFilesIn $fastq_1 $fastq_2 \
  --sjdbGTFfile $gtf \
  --runThreadN $nt \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --quantMode GeneCounts \
  --sjdbGTFtagExonParentGene gene_name \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$sample".sorted"." \
  --outSAMmapqUnique 60 \
  --outSAMattrRGline ID:$RG PU:$platformUnit SM:$sample PL:$platform

# ******************************************
# 1b. Sorting
# No index is produced by star
# bam-file is named back to generic name because it is easier to work with
# ******************************************
module load tools samtools/1.9
mv "$sample".sorted.Aligned.sortedByCoord.out.bam sorted.bam
echo "Loaded samtools 1.9 for indexing"
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
# 2.b. CHECK-POINT FOR NUMBER OF READS
# ******************************************
stats=$(tail -n 2 aln_metrics.txt | head -n 1 | cut -f2,6)
total_reads=$(echo $stats | cut -d ' ' -f1)
pf_reads_aligned=$(echo $stats | cut -d ' ' -f2)
alignment_read_threshold1=1000
alignment_read_threshold2=1000

echo $total_reads were used in alignemnt, $pf_reads_aligned paired reads were aligned ...
if [ $total_reads -lt $alignment_read_threshold1 ]; then
    echo "Not enough reads, exiting ... "
    exit
else
    if [ $pf_reads_aligned -lt $alignment_read_threshold2 ]; then
        echo "Not enough aligned reads, exiting ... "
        exit
        fi
echo "Sufficient paired reads aligned!"
fi
echo "Passed test succesfully!"



# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam

grep "algo: Dedup" -A 4 run.log >> read_summary.txt
reads=$(grep "algo: Dedup" -A 4 run.log | grep 'reads' | cut -d ' ' -f2)
deduplication_threshold=1000

echo $reads reads were used in deduplication ...
if [ $total_reads -lt $deduplication_threshold ]; then
    echo "Not enough reads, exiting ... "
    exit
fi
echo "Passed test succesfully!"


# ******************************************
# 4. Split reads at junction
# This step splits the RNA reads into exon segments by getting rid of Ns while maintaining grouping information, and hard-clips any sequences overhanging into the intron regions.
# Additionally, the step will reassign the mapping qualities from STAR to be consistent with what is expected in subsequent steps by converting from quality 255 to 60.
# second thing is not really necessary (this is from the website / manual), but we keep it to be sure it is done properly
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --algo RNASplitReadsAtJunction --reassign_mapq 255:60 splitted.bam

grep "algo: RNASplitReadsAtJunction" -A 4 run.log >> read_summary.txt
reads=$(grep "algo: RNASplitReadsAtJunction" -A 4 run.log | grep 'reads' | cut -d ' ' -f2)
splitreads_threshold=1000

echo $reads reads left after deduplication ...
if [ $total_reads -lt $deduplication_threshold ]; then
    echo "Not enough reads, exiting ... "
    exit
fi
echo "Passed test succesfully!"



# ******************************************
# 5. Indel realigner
# This step is optional for haplotyper-based caller like HC,
# but necessary for any pile-up based caller.
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i splitted.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels realigned.bam

grep "algo: Realigner" -A 4 run.log >> read_summary.txt
reads=$(grep "algo: Realigner" -A 4 run.log | grep 'reads' | cut -d ' ' -f2)
splitreads_threshold=1000

echo $reads reads after RNASplitReadsAtJunction ...
if [ $total_reads -lt $splitreads_threshold ]; then
    echo "Not enough reads, exiting ... "
    exit
fi
echo "Passed test succesfully!"


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
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper --trim_soft_clip --call_conf 20 --emit_conf 20 -d $dbsnp output.vcf.gz

hc_variant_threshold=-1
nr_variants=$(zgrep -v  '#' output.vcf.gz | wc -l)
echo $nr_variants were called with Haplotyper ...
if [ $total_reads -lt $hc_variant_threshold ]; then
    echo "Not enough variants (less than $hc_variant_threshold) were called!"
    echo "Assuming an error happened, and exiting ... "
    exit
fi
echo "Passed test succesfully!"


# ******************************************
# 6b. Variant calling with TNscope without normal sample
# Here we could include a 'panel of normals' --pon (VCF-file) which should substitute the normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo TNscope --tumor_sample $sample -d $dbsnp --trim_soft_clip output-TNScope.vcf.gz

# ******************************************
# 7. Clean-up and VCF-stats submission
# The final bam-file is the recaled.bam + *bai which may be used in other
# analysis
# VCF-stats are found with vcf statistics and all quality is kept in
# quality_reports
# ******************************************
# rename the final vcf-file to comply with naming scheme
mv output.vcf.gz "$sample".vcf.gz
mv output-TNScope.vcf.gz "$sample"-TNScope.vcf.gz
mv output.vcf.gz.tbi "$sample".vcf.gz.tbi
mv output-TNScope.vcf.gz.tbi "$sample"-TNScope.vcf.gz.tbi
mv run.log "$sample".run.log


$apps/ngs-tools/vcf_statistics.sh "$sample".vcf.gz

# following is just in case quality reports is not named correctly in {vcf/bam}_statistics 
mv quality_reports "$sample".quality_reports

mv *pdf "$sample".quality_reports
mv *txt "$sample".quality_reports
# remove all the files we don't want to keep:
rm recal* 
rm splitted.bam*
rm realigned.bam*
rm score.txt*
rm sorted.bam*
rm deduped.bam*


