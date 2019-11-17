#!/usr/bin/env bash

module load anaconda3/4.4.0
module load bcftools/1.9

# this script collects some stats on the vcf-files AND moves all quality reports to the designated folder 

echo "Using computerome module bcfools/1.9 and python3 modules documented below"

vcf=$1
sample=$(basename $vcf | sed 's/.vcf.*//')
# add optional argument for file naming, default is no prefix 
if [[ -z $2 ]]; then
        opt=$2
        echo "Running with defaults names"
else
    echo Got prefix $2
    opt=$2"_";
fi

start=`date +%s`

echo "# Creating quality_reports dir if it does not already exist"
destination="$sample".quality_reports
mkdir -p $destination

apps="/home/projects/HT2_leukngs/apps/github/code"

# VCF statistics
echo "# Collecting VCF stats for" $1 
bcftools stats $vcf > $destination/"$opt"vcf_summary.txt
grep ^SN -B 1 quality_reports/"$opt"vcf_summary.txt > $destination/"$opt"SN_stats.txt
grep ^QUAL -B 1 quality_reports/"$opt"vcf_summary.txt > $destination/"$opt"QUAL_vcf_stats.txt
grep ^DP -B 1 quality_reports/"$opt"vcf_summary.txt > $destination/"$opt"DP_vcf_stats.txt
echo "# Summarizing and plotting with visualize_stats.py" 
$apps/quality/visualize_stats.py qual $destination/"$opt"QUAL_vcf_stats.txt
$apps/quality/visualize_stats.py dp $destination/"$opt"DP_vcf_stats.txt

# clean up
echo "# Moving files to quality_reports dir" 
mv *metrics* $destination
mv *-report.pdf $destination
mv *summary* $destination

end=`date +%s`
runtime=$((end-start))


echo "# Finished in seconds" $runtime
