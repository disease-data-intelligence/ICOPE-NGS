#!/bin/bash 
#!/bin/bash

input=$1
vcf=$(basename $1 | sed 's/.vep.vcf.gz.*//')
output=$vcf".vep.filter.vcf"
echo -e "Filtering input vcf: $input \t Output name: $output "

module load tools 
module load perl/5.24.0 ensembl-vep/95.1

# Filter on rare alleles (max allele frequency) as we assume negatively selected variants) and variants with high consequences 
# filter_vep --fcrmat vcf - $vcf".vep.vcf.gz" -o $vcf".vep.filter.vcf.gz" --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))"
# only filter on impact and SIFT / PolyPhen 
filter_vep --format vcf -i $input -o $output --force_overwrite --filter "((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))"

nvariants=$(grep -v ^# -c $vcf".vep.filter.vcf")
echo "We get $nvariants variants in the genes of interest"


if [ -f *vep.vcf.gz_summary.html ]; then
	mv *vep.vcf.gz_summary.html quality_reports/vep_summary.html
fi

if [ -f quality_reports/*vep.vcf.gz_summary.html ]; then
	mv quality_reports/*vep.vcf.gz_summary.html quality_reports/vep_summary.html
fi
