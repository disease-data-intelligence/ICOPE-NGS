#!/bin/bash 
#!/bin/bash

vcf=$1 
module load tools 
module load perl/5.24.0 ensembl-vep/95.1

vep -i $vcf".vcf.gz" -o $vcf".vep.vcf.gz" --sift b --polyphen b --symbol --vcf --format vcf --cache --compress_output bgzip --assembly GRCh37 --port 3337 --dir /home/databases/variant_effect_predictor/.vep/

if [ -f *vep.vcf.gz_summary.html ]; then
	mv *vep.vcf.gz_summary.html quality_reports/vep_summary.html
fi

if [ -f quality_reports/*vep.vcf.gz_summary.html ]; then
	mv quality_reports/*vep.vcf.gz_summary.html quality_reports/vep_summary.html
fi
