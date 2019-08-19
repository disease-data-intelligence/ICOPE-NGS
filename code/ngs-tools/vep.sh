#!/bin/bash 
#!/bin/bash
#PBS -W group_list=pr_potentia -A pr_potentia 
#PBS -l nodes=1:ppn=1,mem=50gb,walltime=000:50:00

vcf=$1 
module load tools 
module load perl/5.24.0 ensembl-vep/95.1

vep -i $vcf".vcf.gz" -o $vcf".vep.vcf.gz" --sift b --polyphen b --symbol --vcf --format vcf --cache --compress_output bgzip --assembly GRCh37 --port 3337 --dir /home/databases/variant_effect_predictor/.vep/
