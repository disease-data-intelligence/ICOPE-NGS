#!/bin/bash 
#!/bin/bash

input=$1
vcf=$(basename $1 | sed 's/.vcf.gz.*//')
output=$vcf".vep.vcf.gz"
echo -e "Annotating input vcf: $input \t Output name: $output "

module load tools 
module load perl/5.24.0 ensembl-vep/95.1

vep -i $input -o $output --sift b --polyphen b --symbol --vcf --format vcf --cache --compress_output bgzip --assembly GRCh37 --port 3337 --dir /home/databases/variant_effect_predictor/.vep/

if [ -f *vep.vcf.gz_summary.html ]; then
	mv *vep.vcf.gz_summary.html "$vcf".quality_reports/$vcf"vep_summary.html"
fi

if [ -f quality_reports/*vep.vcf.gz_summary.html ]; then
	mv "$vcf".quality_reports/*vep.vcf.gz_summary.html "$vcf".quality_reports/$vcf"_vep_summary.html"
fi
