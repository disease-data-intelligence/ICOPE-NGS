# Cancer-Seq Exercise 
Cancer-related exercises for the [DTU course 22126](http://teaching.healthtech.dtu.dk/22126/index.php/Program_2020) 

*Adapted from original exercise by Marcin Krzystanek and Aron Eklund.* 
* _Changes include: Some more tools, GATK calls updated to match vs. 4, new resources, comments and questions_    

These exercises will guide you through all steps starting from raw data (FASTQ files)
and resulting in a list of somatic point mutations. 

Estimated time:  2 hours

## Tools 

These exercises are tested with:
* Picard v. 2.21.6 (http://broadinstitute.github.io/picard) 
* GATK  v. 4.1.4.1 (https://github.com/broadinstitute/gatk)
* BWA  v. 0.7.12
* Samtools v. 1.3.1
* Sequenza v. 2.1.2
* R v. 3.5.2 
* TrimGalore v. 0.6.5 
* fastqCombinePairedEnd.py (by [Eric Normandeau](https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py))





## Resources


Known resources are important when working with human data and cancer seq. Luckily, the human is probably one of the most well-annotated species. 
It is generally a good idea to try and use the most up-to-date version of all resources as possible, although a common caveat is that resources have to match and some might not yet be available for the newest versions. 
   
All resources are based on the newest genome build, hg38. The previous and still extensively used build is hg37 (also known as hg19). Matching other resources with the genome build is important as genomic coordinates (location of chromosome e.g. chr:1 differ between genome builds, meaning that for instance a reported SNP has different coordinates depending on the build.

Important points are also that naming needs to match: mentionable naming conventions are defined by ENSEMBL, NCBI (RefSeq) and UCSC. 
In this exercise, the used resources are based on UCSC naming conventions: chromosomes are named chr1, chr2, ..., chrX, chrY and chrM (in "contrast" to ENSEMBL: 1, 2, X, Y and MT, and the less human-readable RefSeq: NC_000001.11, NC_000002.12, NC_000023.11, NC_000024.10 and NC_012920.1). Alternative scaffolds are named chr1_KI270765v1_alt.  
 
 
* GRCh38 (hg38) with UCSC 
* dbSNP 138    
* COSMIC with 30178158 known cancer variants 
* Mills and 1000G genomes Gold Standard indels (known indels)
* gnomAD for allelic frequencies   

All resources are in * . The genome has been bwa-indexed. Consider taking a look at the hg38.fa file. You can look at the header encoding with 

        grep '^>' hg38.fa | head -15    # show first 15   

Besides containing the chromosomes, the genome also contains decoy sequences and HLA-alleles.
## About the data

You will analyze whole-exome sequencing data from a mix of tumor from infiltrating duct adenocarcinoma and head of pancreas, and a matched normal tissue. 
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

The parts where you should actually run the code are marked with a green square: ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 


## PART 1. Raw reads: inspection, QC, cleanup

 

#### 1.1 - Start up ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)  
First, let's define bash variables for keeping our workspace clean.  

        # Input files  
        f1n=/home/27626/exercises/cancer_seq/raw_data/TCRBOA2-N-WEX.read1.fastq.gz
        f2n=/home/27626/exercises/cancer_seq/raw_data/TCRBOA2-N-WEX.read2.fastq.gz
        f1t=/home/27626/exercises/cancer_seq/raw_data/TCRBOA2-T-WEX.read1.fastq.gz
        f2t=/home/27626/exercises/cancer_seq/raw_data/TCRBOA2-T-WEX.read2.fastq.gz
        
        # Resources
        hg38=/home/27626/exercises/cancer_seq/resources/hg38.fa 
        mills=/home/27626/exercises/cancer_seq/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        dbsnp=/home/27626/exercises/cancer_seq/resources/Homo_sapiens_assembly38.dbsnp138.vcf        
        cosmic=/home/27626/exercises/cancer_seq/resources/CosmicCodingMuts_chr_sorted.vcf
        gnomad=home/27626/exercises/cancer_seq/resources/af-only-gnomad.hg38.vcf.gz
        
        # Tools 
        PICARD=/home/27626/bin/picard.jar
        snpSift=/home/27626/bin/snpEff/SnpSift.jar
        # rest are already in your path 

Make a folder for the exercise: 

        mkdir cancer_seq 
        cd cancer_seq 
        outdir=`pwd`
        


#### 1.2 - Take a first look at the data ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 
Look at the data available. You do not need to copy it to your own folder because the files are quite large. 

        ls -l /home/27626/exercises/cancer_seq
        zcat /home/27626/exercises/cancer_seq/TCRBOA2-N-WEX.read1.fastq.gz | head

* How long are the reads? Is your data single or paired end? 
* What type would you prefer for cancer DNA sequencing, and why?
* What is the difference between whole-exome sequencing and RNA-seq? 



#### 1.3 - Read quality trimming and FastQC report (DO NOT RUN)

We do this using [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). 
Trim Galore is actually a wrapper using other trimmers. It can automatically detect adapters for removal, and envoke FastQC simultaneously.  


        ### Trim reads with trim_galore wrapper, produce both fastqc and trimming reports
        trim_galore --fastqc --fastqc_args "--outdir ${outdir}/trimmed_normal" --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1n $f2n
        trim_galore --fastqc --fastqc_args "--outdir ${outdir}/trimmed_tumor" --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1t $f2t


* What does the argument `--quality 20` mean? Get help by running:
        
        trim_galore --help


#### 1.4 - Keeping files in sync (DO NOT RUN)
Trim Galore does not make sure that all reads are paired which is a requirement for bwa. 
This is handled with fastqCombinePairedEnd.py. 
 
    
        python fastqCombinePairedEnd.py $outdir/trimmed_normal/TCRBOA2-N-WEX.read1_trimmed.fq.gz $outdir/trimmed_normal/TCRBOA2-N-WEX.read2_trimmed.fq.gz
        python fastqCombinePairedEnd.py $outdir/trimmed_normal/TCRBOA2-N-WEX.read1_trimmed.fq.gz $outdir/trimmed_normal/TCRBOA2-N-WEX.read2_trimmed.fq.gz



Set up new variables for the newly created files. 

        f1n_val=$outdir/trimmed_normal/TCRBOA2-N-WEX.read1_trimmed.fq.gz_pairs_R1.fastq.gz
        f2n_val=$outdir/trimmed_normal/TCRBOA2-N-WEX.read2_trimmed.fq.gz_pairs_R2.fastq.gz
        f1t_val=$outdir/trimmed_tumor/TCRBOA2-T-WEX.read1_trimmed.fq.gz_pairs_R1.fastq.gz
        f2t_val=$outdir/trimmed_tumor/TCRBOA2-T-WEX.read2_trimmed.fq.gz_pairs_R2.fastq.gz


## PART 2. Alignment and additional preprocessing (DO NOT RUN)

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

        # make a directory for alignment and enter 
        mkdir aligned
        cd aligned 
           
        ### Run bwa mem
        bwa mem -M -t 4 -R $ReadGoupID_N $hg38 $f1n_val $f2n_val \
            | samtools view -Sb -@ 1 - | samtools sort -@ 3 > TCRBOA2-N-WEX.bam 
        bwa mem -M -t 4 -R $ReadGoupID_T $hg38 $f2t_val $f2t_val \
            | samtools view -Sb -@ 1 - | samtools sort -@ 3 > TCRBOA2-T-WEX.bam

Optionally, the @RG line can provide additional information; please see the 
[SAM format specification](http://www.samformat.info) as well as [samtools webpage](http://samtools.sourceforge.net) if you want to know more.
 
The command after the pipe is to compress the output of bwa (sam-format) to binary (bam-format).  
We also sorting the bam-file to avoid unnecessary intermediate files, keeping in mind a bam-file might take up a lot of space. 
This can also be done separately with the following command: 
  
Note that it is also possible to pipe the output of bwa directly to `samtools` to sort (without writing an unsorted bam-file).

        samtools sort -@ 3 TCRBOA2-N-WEX_unsorted.bam -o TCRBOA2-N-WEX.bam
        samtools sort -@ 3 TCRBOA2-T-WEX_unsorted.bam -o TCRBOA2-N-WEX_unsorted.bam


#### 2.2 - Mark duplicates (DO NOT RUN) 
_~40/~58 min per file_ 

We use [Picard](https://broadinstitute.github.io/picard/) to mark PCR duplicates so that
they will not introduce false positives and bias in the subsequent analysis.

        mkdir tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=TCRBOA2-N-WEX.bam OUTPUT=TCRBOA2-N-WEX_deduped.bam METRICS_FILE=TCRBOA2-N-WEX.metrics.txt \
            TMP_DIR=./tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=TCRBOA2-T-WEX.bam OUTPUT=TCRBOA2-T-WEX_deduped.bam METRICS_FILE=TCRBOA2-T-WEX.metrics.txt \
            TMP_DIR=./tmp        
        
      

#### 2.3 - Index the BAM files (DO NOT RUN) (3.5 min pr. file)
_~3 min per file_

        samtools index TCRBOA2-T-WEX_deduped.bam
        samtools index TCRBOA2-N-WEX_deduped.bam



#### 2.4a - BaseRecalibrator - Part 1 (DO NOT RUN)
_56/53 min per file_

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
        
        
#### 2.4b - Apply BaseRecalibration - Part 2 (DO NOT RUN)
_~34 and ~32 min per file_ 

         gatk ApplyBQSR \
           -I TCRBOA2-N-WEX_deduped.bam \
           -R $hg38 \
           --bqsr-recal-file normal.recal.table \
           -O TCRBOA2-N-WEX_recaled.bam
          
         gatk ApplyBQSR \
           -I TCRBOA2-T-WEX_deduped.bam \
           -R $hg38 \
           --bqsr-recal-file tumor.recal.table \
           -O TCRBOA2-T-WEX_recaled.bam
        

Now, the resulting BAM files are ready to be processed with MuTect2.


## PART 3. Somatic mutation calling (BAM file -> VCF file)
 

#### 3.1 - MuTect2 ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 
_All chromosomes took ~327 minutes, (~5.4 hours with four threads)_

We are going to use "soft-links" to the input bam-files. This a way to avoid working with long paths and is an alternative to defining bash-variables. 

        mkdir variant_calling 
        cd variant_calling
        
        # create links   
        ln -s /home/27626/exercises/cancer_seq/aligned/TCRBOA2-T-WEX_recaled.bam . 
        ln -s /home/27626/exercises/cancer_seq/aligned/TCRBOA2-N-WEX_recaled.bam . 
       
Use `ls -l` to see how they appear in your directory. 


We use [MuTect2][MuTect2], a somatic mutation caller that identifies both SNV and indels. 
The produced VCF-file (variant calling format) is the standard format for storing variants, although the output of Mutect2 has some information specific for somatic variants. See here for specs: https://samtools.github.io/hts-specs/VCFv4.1.pdf
   
A big difference in cancer-seq variant calling using Mutect2 is that there are no ploidy assumptions. 
This accommodates tumor data which can have many copy number variants (CNVs).   


[MuTect2]: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php

Mutect2 is computationally intensive so we recommend to parallelize if possible. 
One way to achieve this is to split processes by chromosomes (calling variants for each chromosome and then merging vcf-files.) 

Since we do not have the time and capacity to process the entire genome during our
exercises, we will call somatic mutations on a small part of chromosome 10, 
from the 3,100,000th to the 5,100,000th base pair, which is set with the `-L` option. 

Now we are ready for finding variants! 


        # Set chromosome and location:
        CHR_LOC=chr10:3100000-5100000
        # Run Mutect2
        gatk Mutect2 \
            -R $hg38 \
            -I TCRBOA2-N-WEX_recaled.bam \
            -I TCRBOA2-T-WEX_recaled.bam \
            -normal TCRBOA2-N-WEX \
            -L ${CHR_LOC} \
            --germline-resource $gnomad \
            -O TCRBOA2_${CHR_LOC}.vcf
            
        ### To process the whole genome, simply omit the -L option.

Take a look at the resulting VCF file with `less -RS TCRBOA2_${CHR_LOC}.vcf` Then try to count the number of raw variants: 

         grep -v "^#" TCRBOA2_${CHR_LOC}.vcf | wc -l

* How many variants did you find?   


#### 3.2 Compare with calling all variants ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)  
Just for comparison, we try to call all variants in the interval for the germline and the tumor sample.

        gatk HaplotypeCaller -I TCRBOA2-T-WEX_recaled.bam -R $hg38 -L ${CHR_LOC} -O TCRBOA2-T.vcf --dbsnp $dbsnp 
        gatk HaplotypeCaller -I TCRBOA2-N-WEX_recaled.bam -R $hg38 -L ${CHR_LOC} -O TCRBOA2-N.vcf --dbsnp $dbsnp 

Count the variant lines in each with the same command as above: 

        grep -v "^#" TCRBOA2-T.vcf | wc -l
        grep -v "^#" TCRBOA2-N.vcf | wc -l 
        

* Where do you find the most raw variants and does that make biological sense? What is the difference in the two numbers and does it match above?  

#### 3.3 - Filter the VCF output ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 
Before continuing, we need to filter the raw vcf-output to only get confident variants:  
 
       gatk FilterMutectCalls \
          -V TCRBOA2_${CHR_LOC}.vcf \
          -R $hg38 \
          -O TCRBOA2_${CHR_LOC}_filtered.vcf 

Try to look at the output with `less -RS TCRBOA2_${CHR_LOC}_filtered.vcf`. 
* What does it look like in the filter column? What kind of filters were applied? 


To add some extra information to the vcf-file, we will also annotate with SNP-ids. HaplotypeCaller can do this as it calls variants, but using Mutect2 we need to do it ourselves: 

        java -jar $snpSift annotate $dbsnp TCRBOA2_${CHR_LOC}_filtered.vcf > TCRBOA2_${CHR_LOC}_filtered_anno.vcf 

 
Now try to filter mutational calls by selecting those with Mutect "PASS" annotation.

        grep PASS TCRBOA2_${CHR_LOC}_filtered_anno.vcf | grep -v "^#"  

You should at least see this line (without the header). Don't forget you can scroll to the sides! 

        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCRBOA2-N-WEX   TCRBOA2-T-WEX
        chr10	3165513	rs9423502	G	C	.	PASS	CONTQ=93;DP=168;ECNT=1;GERMQ=93;MBQ=38,39;MFRL=225,184;MMQ=60,60;MPOS=4;NALOD=1.86;NLOD=21.37;POPAF=1.32;SEQQ=19;STRANDQ=16;TLOD=6.08;CAF=[0.9454,0.05464];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;NSM;OTHERKG;PH3;REF;RS=9423502;RSPOS=3207705;S3D;SAO=0;SLO;SSR=0;VC=SNV;VLD;VP=0x050300000a01150517000101;WGT=1;dbSNPBuildID=119	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:79,0:0.014:79:41,0:38,0:19,60,0,0	0/1:79,3:0.054:82:41,2:38,1:20,59,0,3

There is a lot of information in the line so no surprise of you're feeling a bit overwhelmed. A brief explanation of each part of the line above is in the header of the VCF file 
(use `less -RS TCRBOA2_${CHR_LOC}_filtered.vcf` to look at it, but it is also difficult to interpret). 

Let's focus on the FORMAT lines. The information kept here is organized by `GT:AD:AF:DP:F1R2:F2R1:SB` and the values are kept in the two proceeding columns, also seperated by colons. 
The first starting with 0/0 refers to the normal sample, whereas the column beginning with 0/1 refers to the tumor. 
This means that the tumor is heterozygote (REF/ALT) for the mutation which is not seen in the germline sample at all (it is REF/REF).  

After genotype (GT) we have allelic depth (AD) which is "79,3" (i.e. 79 reads and 3 reads for
the reference and mutant allele respectively) in the tumor and "79,0" in the normal sample. Then comes allelic frequency, which is
a fraction of the mutant allele out of all aligned bases in this position and the depth. We will skip the remaining values for F1R2, F2R1 and SB for now. 



## PART 4. Interpretation of the resulting somatic mutations  

A list of chromosome coordinates is kind of hard to interpret. Here are some ways to approach the results.  

#### 4.1 - Variant annotation with dbSNP, Variant effect predictor and COSMIC  ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 


* Find the RS identifier from the cancer mutation and look it up at [dbSNP](https://www.ncbi.nlm.nih.gov/snp). 
    * Find the frequency table tab. Is your mutation common in some populations? What does a high frequency tell you about its role in cancer? 



So far you have processed and analyzed only a small section of chromosome 10.

Now, let's analyze a bigger piece of the genome. Pick your favorite chromosome and
find the corresponding VCF file on the server.  For example, if you choose
chromosome 7, you would use this file:

        ls -l /home/27626/exercises/cancer_seq/variant_calling/split_by_chromosome/chr7.vcf

Hint: your results will be more interesting if you pick chromosome 
6, 13, 15, 17, 19, 20, 22, or X!

Filter the VCF to retain only the lines marked as "PASS".  

        grep "PASS" /home/27626/exercises/cancer/chr7.vcf > filtered.chr7.vcf

* Download the *filtered* VCF to your own computer using `scp ngsstudXXX@10.44.152.8:/path/to/file . ` and submit it to the 
[VEP website](http://www.ensembl.org/Tools/VEP) using default settings. 
When the results become available, look in the "Somatic status" column. Are there
any known cancer mutations? 

* If you find a known cancer mutation, find its COSMIC identifier 
(COSM######, e.g. COSM4597270) in the "existing variant" column.
Search for your COSMIC identifier in the
[COSMIC database](http://cancer.sanger.ac.uk/cosmic).
    * In which tissues is this mutation found?

* Go back to the VEP output and find the columns SIFT and PolyPhen (consequence predictors). Sort your results by either of the scores and see if you find any damaging or deleterious scores. 
    * What is the difference between how SIFT and PolyPhen predicts consequences? [Look here for a hint](https://m.ensembl.org/info/genome/variation/prediction/protein_function.html)
    * Do SIFT and PolyPhen agree? 


#### 4.2 - cBioPortal  ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 

Go to [cBioPortal](http://www.cbioportal.org), a website that provides tools to analyze
several large cancer sequencing datasets. In Quick Select, choose "TCGA PanCancer Atlas studies".
Then press "Query by Gene" and type in the name of the gene that was hit by this
mutation. Choose "mutations" as we have not looked at Copy Number Alterations. Press "Submit Query". 
Look at the barcharts and play around with the options.
* How often is this gene mutated in various cancer types?  



#### 4.3 - Inference of tissue of origin ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) 

Next we'll do some analysis on a VCF file containing somatic mutations found throughout
the entire genome:

        ls /home/27626/exercises/cancer_seq/variant_calling/TCRBOA2.vcf 

Unlike VEP, TumorTracer requires VCF files to have the header information.
Thus, we will filter this VCF file to retain: 1) header lines (which begin with "#"),
and 2) data lines with a PASS call.
        
        grep -E "^#|PASS" /home/27626/exercises/cancer_seq/variant_calling/TCRBOA2.vcf > TCRBOA2_filtered.vcf

* Download the vcf-file and submit it to the 
[TumorTracer server](http://www.cbs.dtu.dk/services/TumorTracer/).
Make sure to specify that this VCF was generated using GRCh38 coordinates.

* What tissue does TumorTracer predict?  Is it a confident prediction?

