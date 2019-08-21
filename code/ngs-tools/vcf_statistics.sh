#!/usr/bin/env bash

module load bcftools/1.9

# this script collects some stats on the vcf-files AND moves all quality reports to the designated folder 

echo "Using computerome module bcfools/1.9 and python3 modules documented below"

vcf=$1
start=`date +%s`

echo "# Creating quality_reports dir if it does not already exist" 
mkdir -p quality_reports 

apps="/home/projects/HT2_leukngs/apps/github/code"

# VCF statistics
echo "# Collecting VCF stats"
bcftools stats $vcf > quality_reports/vcf_summary.txt
grep ^SN -B 1 quality_reports/vcf_summary.txt > quality_reports/SN_stats.txt
grep ^QUAL -B 1 quality_reports/vcf_summary.txt > quality_reports/QUAL_vcf_stats.txt
grep ^DP -B 1 quality_reports/vcf_summary.txt > quality_reports/DP_vcf_stats.txt
echo "# Summarizing and plotting with visualize_stats.py" 
$apps/quality/visualize_stats.py qual quality_reports/QUAL_vcf_stats.txt
$apps/quality/visualize_stats.py dp quality_reports/DP_vcf_stats.txt

# clean up
echo "# Moving files to quality_reports dir" 
mv *metrics* quality_reports
mv *-report.pdf quality_reports
mv *summary* quality_reports

end=`date +%s`
runtime=$((end-start))


echo "# Finished in seconds" $runtime
