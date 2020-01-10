# Cancer-Seq Exercise 
Cancer-related exercises for the [DTU course 22126](http://teaching.healthtech.dtu.dk/22126/index.php/Program_2020) 

*Adapted from original exercise by Marcin Krzystanek and Aron Eklund.* 

These exercises will guide you through all steps starting from raw data (FASTQ files)
and resulting in a list of somatic point mutations. 

Estimated time:  2 hours

## Prerequisites

These exercises are tested with:
* Picard v. 2.21.6 “Picard Toolkit.” 2019. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/; Broad Institute 
* GATK  v. 4.1.4.1 (https://github.com/broadinstitute/gatk)
* BWA  v. 0.7.12
* Samtools v. 1.3.1
* Sequenza v. 2.1.2
* R v. 3.5.2 
* TrimGalore v. 0.6.5 
* fastqCombinePairedEnd.py by [Eric Normandeau](https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py)


https://samtools.github.io/hts-specs/VCFv4.1.pdf


## Resources


Known resources are important when working with human data and cancer seq. Luckily, the human is probably one of the most well-annotated species. 
It is generally a good idea to try and use the most up-to-date version of all resources as possible, although a common caveat is that resources have to match and some might not yet be available for the newest versions. 
   
All resources are based on the newest genome build, hg38. The previous and still extensively used build is hg37 (also known as hg19). Matching other resources with the genome build is important as genomic coordinates (location of chromosome e.g. chr:1 differ between genome builds, meaning that for instance a reported SNP has different coordinates depending on the build.

Important points are also that naming needs to match: mentionable naming conventions are defined by ENSEMBL, NCBI (RefSeq) and UCSC. 
In this exercise, the used resources are based on UCSC naming conventions: chromosomes are named chr1, chr2, ..., chrX, chrY and chrM (in "contrast" to ENSEMBL: 1, 2, X, Y and MT, and the less human-readable RefSeq: NC_000001.11, NC_000002.12, NC_000023.11, NC_000024.10 and NC_012920.1). Alternative scaffolds are named chr1_KI270765v1_alt.  
 
 
* GRCh38 (hg38) with UCSC 
* dbSNP  
* COSMIC 30178158 known cancer variants 
* Mills and 1000G genomes Gold Standard indels (known indels)  

All resources are in * . The genome has been bwa-indexed. Consider taking a look at the hg38.fa file. You can look at the header encoding with 

        grep '^>' hg38.fa | head -15    # show first 15   

Besides containing the chromosomes, the genome also contains decoy sequences, HLA-alleles and 
## About the data

You will analyze whole-exome sequencing data from a mix of infiltrating duct adenocarcinoma and head of pancreas, and a matched normal tissue. 
This is known as somatic variant calling or as paired variant calling because we have a tumor-normal pair, and germline variation will be filtered from our final set of variants. 
Thus we will only find mutations that are specific for the tumor and potential _driver mutations_.  

The data used in this exercise has been released for scientific and educational use by the
[Texas Cancer Research Biobank](http://txcrb.org/data.html) and is fully described in 
[this paper](https://www.nature.com/articles/sdata201610). You can read about the sequencing protocol [here](http://txcrb.org/sequencing.html). 


Please note the Conditions of Data Use:

> By downloading or utilizing any part of this dataset, end users must agree to the following conditions of use:
> * No attempt to identify any specific individual represented by these data or any derivatives of these data will be made.
> * No attempt will be made to compare and/or link this public data set or derivatives in part or in whole to private health information.
> * These data in part or in whole may be freely downloaded, used in analyses and repackaged in databases.
> * Redistribution of any part of these data or any material derived from the data will include a copy of this notice.
> * The data are intended for use as learning and/or research tools only.
> * This data set is not intended for direct profit of anyone who receives it and may not be resold.
> * Users are free to use the data in scientific publications if the providers of the data (Texas Cancer Research Biobank and Baylor College of Medicine Human Genome Sequencing Center) are properly acknowledged.

The raw data files are located on the server at `/home/27626/exercises/cancer_seq`

Much of this exercise is based on the Best Practices by [Broad Institute](https://www.broadinstitute.org/) by  MIT and Harvard and located in Cambridge, Massachusetts.
The Broad provides the Genome Analysis Toolkit (gatk) which you will probably become very familiar with if you continue with genomics and NGS in the future. 
Best practices for somatic variant calling can be found [here](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146).
   
  
## Somatic point mutation exercise

**IMPORTANT IMPORTANT IMPORTANT** - Since the full procedure takes a long time, 
we will **not** ask you to perform the full alignment and full mutation calling, and assume that you have become familiar with this in previous exercises.  
However, for reference, we provide the code needed for the full analysis.  
Thus, you can use this code later in the course project or in your own work, 
should you work with cancer patient sequencing data.

The parts where you should actually run the code include: 1.1, 1.2, 3, and 4


### PART 1. Raw reads: inspection, QC, cleanup

#### 1.1 - Take a first look at the data

        ls /home/27626/exercises/cancer_seq
        zcat /home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read1.fastq.gz | head

Q1: How long are the reads? Is your data single or paired end? 
What type would you prefer for cancer DNA sequencing, and why?
What is the difference between whole-exome sequencing and RNA-seq? 


#### 1.2 - Define some bash variables

        ### Define bash variables for brevity:
        f1n=/home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read1.fastq.gz
        f2n=/home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read2.fastq.gz
        f1t=/home/27626/exercises/cancer_seq/TCRBOA2-T-WEX.read1.fastq.gz
        f2t=/home/27626/exercises/cancer_seq/TCRBOA2-T-WEX.read2.fastq.gz
        HREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set
        FREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fa 
        mills=/home/27626/exercises/cancer/human_GRCh38/Indel_refs/mills_gold.b38.vcf
        SREFF=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/1000G.snps.b38.vcf
        dbsnp=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/dbsnp/All_20160527_chr.vcf
        cosmicREFF=/home/27626/exercises/cancer/human_GRCh38/cosmic/CosmicCodingMuts_chr_sorted.vcf
        GATK=/home/27626/exercises/cancer/programs/GenomeAnalysisTK.jar
        PICARD=/home/27626/bin/picard.jar
        # SAMTOOLS=/home/27626/bin/samtools # this works wihtout path 
        TRIM_GALORE=/home/27626/exercises/cancer/programs/trim_galore
        outdir=`pwd`


#### 1.3 - Read quality trimming and FastQC report (DO NOT RUN)

We do this using [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). 
Trim Galore is much like other trimmers and can automatically detect adapters for removal, and envoke FastQC simultaneously.  


        ### Arguments to be passed to FastQC
        args="'--outdir ${outdir}'"
        
        ### Trim reads with trim_galore wrapper, produce both fastqc and trimming reports
        trim_galore --fastqc --fastqc_args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired /home $f2n
        trim_galore --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1t $f2t

        trim_galore --fastqc --fastqc_args "--outdir test" --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired /home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read1.fastq.gz /home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read2.fastq.gz 


Q2: What does the argument `--quality 20` mean? Get help by running:
        
        trim_galore --help



#### 1.4 - Keeping files in sync (DO NOT RUN)
Trim Galore does not make sure that all reads are paired which is a requirement for bwa. 
This is handled with fastqCombinePairedEnd.py. 
 
 
    python fastqCombinePairedEnd.py ~/texas/trimmed_normal/TCRBOA2-N-WEX.read1_trimmed.fq.gz ~/texas/trimmed_normal/TCRBOA2-N-WEX.read2_trimmed.fq.gz
    python fastqCombinePairedEnd.py ~/texas/trimmed_normal/TCRBOA2-N-WEX.read1_trimmed.fq.gz ~/texas/trimmed_normal/TCRBOA2-N-WEX.read2_trimmed.fq.gz



Set up new variables for the newly created files. I assume the validated and filtered
files were created in your working directory (for me this is /home/27626/exercises/cancer/
so you can find these files there if you need them).

        f1n_val=TCRBOA2-N-WEX.read1.fastq.bz2_val_1.fq.gz
        f2n_val=TCRBOA2-N-WEX.read2.fastq.bz2_val_2.fq.gz
        f1t_val=TCRBOA2-T-WEX.read1.fastq.bz2_val_1.fq.gz
        f2t_val=TCRBOA2-T-WEX.read2.fastq.bz2_val_2.fq.gz


### PART 2. Alignment and additional preprocessing (DO NOT RUN)

#### 2.1 - Alignment (DO NOT RUN) 
_~4.5 hours with 4 processors, ~37 minutes with 14 processors_

We use [bwa mem](https://github.com/lh3/bwa) for aligning reads to the genome. 
We align the tumor sample and normal sample separately. 

Importantly, a Read Group ID line (@RG line) must be defined by the user, because Mutect2
and other programs in the pipeline below depend on information in this line. Here we
demonstrate one way of adding the @RG line to the resulting BAM file:

        ### @RG ID # read group ID, needs to be unique for fastq file due to downstream processing, takes\
        preference when used by some programs
        ### @RG SM # sample ID, unique for each tumor and normal sample, not to be confused with patient ID
        ### @RG PL # platform name
        ### @RG LB # library name
        ### @RG PU # Platform unit, needs to be unique for fastq file due to downstream processing, takes\
        preference when used by some programs
        ### Let's create an @RG line that we will use when running bwa mem alignment
        ReadGoupID_N="\"@RG\tID:TCRBOA2-N-WEX\tSM:TCRBOA2-N-WEX\tPL:ILLUMINA\tLB:libN\tPU:TCRBOA2-N-WEX"\"
        ReadGoupID_T="\"@RG\tID:TCRBOA2-T-WEX\tSM:TCRBOA2-T-WEX\tPL:ILLUMINA\tLB:libT\tPU:TCRBOA2-T-WEX"\"

        ### Run bwa mem
        bwa mem -M -t 4 -R $ReadGoupID_N $HREFF $f1n_val $f2n_val \
            | samtools view -Sb -@ 1 - | samtools sort -@ 3 > patient2_n.bam  
        bwa mem -M -t 4 -R $ReadGoupID_T $HREFF $f2t_val $f2t_val \
            | samtools view -Sb -@ 1 - | samtools sort -@ 3 > patient2_t.bam

Optionally, the @RG line can provide additional information; please see the 
[SAM format specification](http://www.samformat.info) as well as [samtools webpage](http://samtools.sourceforge.net) if you want to know more.
 
The command after the pipe is to compress the output of bwa (sam-format) to binary (bam-format).  
We also sorting the bam-file to avoid unnecessary intermediate files, keeping in mind a bam-file might take up a lot of space. 
This can also be done separately with the following command: 
  
Note that it is also possible to pipe the output of bwa directly to samtools to sort (without writing an unsorted bam-file).

        samtools sort -@ 3 patient2_n.bam -o patient2_n.sorted.bam
        samtools sort -@ 3 patient2_t.bam -o patient2_t.sorted.bam


#### 2.3 - Mark duplicates (DO NOT RUN) 
_~40/~58 min pr file_ 

We use [Picard](https://broadinstitute.github.io/picard/) to mark PCR duplicates so that
they will not introduce false positives and bias in the subsequent analysis.

        mkdir tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=patient2_n.sorted.bam OUTPUT=patient2_n.sorted.dedup.bam METRICS_FILE=patient2_n.metrics.txt \
            TMP_DIR=./tmp
        java -Xmx10G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar $PICARD MarkDuplicates\
            INPUT=patient2_t.sorted.bam OUTPUT=patient2_t.sorted.dedup.bam METRICS_FILE=patient2_t.metrics.txt \
            TMP_DIR=./tmp
        
        
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar ../../picard/build/libs/picard.jar MarkDuplicates INPUT=TCRBOA2-N-WEX.bam OUTPUT=TCRBOA2-N-WEX_deduped.bam METRICS_FILE=patient2_n.metrics.txt

#### 2.4 - Index the BAM files (DO NOT RUN) (3.5 min pr. file)
_3 min per file_

        samtools index patient2_n.sorted.dedup.bam
        samtools index patient2_t.sorted.dedup.bam



#### 2.5 - BaseRecalibrator - Part 1 (DO NOT RUN)
_56 & 53 min pr file_

We use [GATK](https://software.broadinstitute.org/gatk/) to recalibrate base quality scores. 

Each base in each sequencing read comes out of the sequencer with an individual quality score. 
Depending on the machine used for sequencing, these scores are subject to various
sources of systematic technical error. Base quality score recalibration (BQSR) works by
applying machine learning to model these errors empirically and adjust the quality scores
accordingly. 
Read more about it [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)

Link to [tool documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036712791-BaseRecalibrator).
--known-sites specifies one or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.

There is more information on BSQR [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php). 

        gatk BaseRecalibrator \
            -I TCRBOA2-N-WEX_deduped.bam \
            -R $hg38 \
            --known-sites $dbsnp \
            --known-sites $mills \
            -O normal.recal.table
            
        gatk BaseRecalibrator \
            -I TCRBOA2-T-WEX_deduped.bam \
            -R $hg38 \
            --known-sites $dbsnp \
            --known-sites $mills \
            -O tumor.recal.table          
        
        
#### 2.6 - Apply BaseRecalibration - Part 2 (DO NOT RUN)
_~34 and ~32 min per file_ 

         gatk ApplyBQSR \
           -I TCRBOA2-N-WEX_deduped.bam \
           -R $hg38 \
           --bqsr-recal-file normal.recal.table \
           -O TCRBOA2-N-WEX_final.bam
          
         gatk ApplyBQSR \
               -I TCRBOA2-T-WEX_deduped.bam \
               -R $hg38 \
               --bqsr-recal-file tumor.recal.table \
               -O TCRBOA2-T-WEX_final.bam
        

Now, the resulting BAM files are ready to be processed with MuTect2.


### PART 3. Somatic mutation calling (BAM file -> VCF file)
_All sequences: ~327 minutes, (~5.4 hours_ 

#### 3.1 - MuTect2
# https://github.com/broadinstitute/gatk/issues/4390
We use [MuTect2][MuTect2], a somatic mutation caller that identifies both SNV and indels. 
The produced VCF-file (variant calling format) is the standard format for storing variants, although the output of Mutect2 has some information specific for somatic variants.   
A big difference in cancer-seq variant calling using Mutect2 is that, although it uses the regular variant caller, HaplotypeCaller’s approach, there are no ploidy assumptions. 
This accommodates tumor data which can have many copy number variants (CNVs).   


[MuTect2]: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php

Mutect2 is computationally intensive so we recommend to parallelize if possible. 
One way to achieve this is to split processes by chromosomes (calling variants for each chromosome and then merging vcf-files.) 

Since we do not have the time and capacity to process the entire genome during our
exercises, we will call somatic mutations on a small part of chromosome 1, 
from the 50.000.000th to the 52.000.000th base pair.

        # Set chromosome and location:
        CHR_LOC=chr1:50000000-52000000
        # Run Mutect2
        gatk Mutect2 \
            -R ../../references/hg38.fa \
            -I TCRBOA2-N-WEX_recaled.bam \
            -I TCRBOA2-T-WEX_recaled.bam \
            -normal TCRBOA2-N-WEX \
            -L ${CHR_LOC} \
            --germline-resource $gnomad \ 
            -O TCRBOA2_somatic_${CHR_LOC}.vcf
            
        ### To process the whole genome, simply omit the -L option.

Take a look at the resulting VCF file.  


#### 3.2 - Filter the VCF output
Just like with HaplotypeCaller, we need to filter the raw vcf-output. 
 
       gatk FilterMutectCalls \
          -V TCRBOA2.vcf \
          -R $hg38 \
          -O TCRBOA2_filtered.vcf 


For a start, try to filter mutational calls by selecting those with MuTect "PASS" annotation.

        cat patient2_t.${CHR_LOC}.mutect2.vcf | grep PASS

You should see this line:
"chr1	50973993	rs746646631	C	T	.	PASS	DB;ECNT=1;HCNT=2;MAX_ED=.;MIN_ED=.;NLOD=33.99;TLOD=7.29	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1	0/1:129,6:0.044:3:3:0.500:3973,169:66:63	0/0:132,0:0.00:0:0:.:4093,0:75:57"

A brief explanation of each part of the line above is in the header of the VCF file 
(use "less -RS patient2_t.${CHR_LOC}.mutect2.vcf" to look at it). Importantly, the column starting with 0/0 refers to the
normal sample, whereas the column beginning with 0/1 refers to the tumor. 
After genotype (GT) we have allelic depth (AD) which is "129,6" (i.e. 129 and 6 for
the reference and mutant allele respectively). Then comes allelic frequency, which is
a fraction of the mutant allele out of all aligned bases in this position. 
For more information about the MuTect2 output go to [MuTect2][MuTect2].


### PART 4. Interpretation of the resulting somatic mutations

#### 4.1 - Interpretation of somatic mutations

Q3: Try to search [dbSNP](https://www.ncbi.nlm.nih.gov/snp) for rs746646631. 
What gene does it belong to? Is this mutation protein-changing?

Go to [cBioPortal](http://www.cbioportal.org), a website that provides tools to analyze
several large cancer sequencing datasets. Type the name of the gene that was hit by this
mutation in the "Enter Gene Set:" box in the bottom of the page and press submit. 
How often is this gene mutated in various cancer types?  


#### 4.2 - More interpretation of somatic mutations

So far you have processed and analyzed only a small section of chromosome 1.

Now, let's analyze a bigger piece of the genome. Pick your favorite chromosome and
find the corresponding VCF file on the server.  For example, if you choose
chromosome 7, you would use this file:

        ls /home/27626/exercises/cancer/patient2.chr7.mutect2.vcf

Hint: your results will be more interesting if you pick chromosome 
6, 13, 15, 17, 19, 20, 22, or X!

Filter the VCF to retain only the lines marked as "PASS".  

        grep "PASS" /home/27626/exercises/cancer/patient2.chr7.mutect2.vcf > filtered.chr7.vcf

Download the *filtered* VCF to your own computer and submit it to the 
[VEP website](http://www.ensembl.org/Tools/VEP) using default settings. 
When the results become available, look in the "Somatic status" column. Are there
any known cancer mutations?
If you find a known cancer mutation, find its COSMIC identifier 
(COSM######, e.g. COSM4597270) in the "existing variant" column.
Search for your COSMIC identifier in the
[COSMIC database](http://cancer.sanger.ac.uk/cosmic).
In which tissues is this mutation found?


#### 4.3 - Inference of tissue of origin

Next we'll do some analysis on a VCF file containing somatic mutations found throughout
the entire genome:

        ls /home/27626/exercises/cancer/patient2_t.mutect2.vcf

Unlike VEP, TumorTracer requires VCF files to have the header information.
Thus, we will filter this VCF file to retain: 1) header lines (which begin with "#"),
and 2) data lines with a PASS call.
        
        grep -E "^#|PASS" /home/27626/exercises/cancer/patient2_t.mutect2.vcf > filtered.patient2_t.mutect2.vcf

Submit the filtered VCF to the
[TumorTracer server](http://www.cbs.dtu.dk/services/TumorTracer/).
Make sure to specify that this VCF was generated using GRCh38 coordinates.

What tissue does TumorTracer predict?  Is it a confident prediction?

