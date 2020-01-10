# Notes and Resources 

## Basics 

Download stuff on linux/Computerome with  

        wget /path/to/file       # just get file 
        wget -r /path/to/file   # copy the folder structure of source  


### Naming stuff 
Important points are also that naming needs to match: mentionable naming conventions are defined by ENSEMBL, NCBI (RefSeq) and UCSC. 
        
        # Assembly name:  GRCh38.p12
        # Description:    Genome Reference Consortium Human Build 38 patch release 12 (GRCh38.p12)
        # Organism name:  Homo sapiens (human)
        # Taxid:          9606
        # BioProject:     PRJNA31257
        # Submitter:      Genome Reference Consortium
        # Date:           2017-12-21
        # Assembly type:  haploid-with-alt-loci
        # Release type:   patch
        # Assembly level: Chromosome
        # Genome representation: full
        # RefSeq category: Reference Genome
        # GenBank assembly accession: GCA_000001405.27
        # RefSeq assembly accession: GCF_000001405.38


Some explanations and discussions around this [here](https://www.biostars.org/p/210935/), [here](https://www.biostars.org/p/314017/) and [here](https://www.biostars.org/p/396813/).  
GRCh38.p12 means that it is a patch of GRCh38 and should have matching annotations. Read about that [here](https://www.biostars.org/p/376195/). 

Example of conversions in conversion table found here::  

        ## GenBank Unit Accession       RefSeq Unit Accession   Assembly-Unit name
        GCA_000001305.2      GCF_000001305.15        Primary Assembly

        
 
### Example of name conversions
Chromosomes are named chr1, chr2, ..., chrX, chrY and chrM. 
 
        
        # Sequence-Name Sequence-Role   Assigned-Molecule       Assigned-Molecule-Location/Type GenBank-Accn    Relationship    RefSeq-Accn     Assembly-Unit   Sequence-Length UCSC-style-name
        1       assembled-molecule      1       Chromosome      CM000663.2      =       NC_000001.11    Primary Assembly        248956422       chr1
        X       assembled-molecule      X       Chromosome      CM000685.2      =       NC_000023.11    Primary Assembly        156040895       chrX
        Y       assembled-molecule      Y       Chromosome      CM000686.2      =       NC_000024.10    Primary Assembly        57227415        chrY
        MT      assembled-molecule      MT      Mitochondrion   J01415.2        =       NC_012920.1     non-nuclear     16569   chrM
        HSCHR3UN_CTG2   unlocalized-scaffold    3       Chromosome      GL000221.1      =       NT_167215.1     Primary Assembly        155397  chr3_GL000221v1_random
        HSCHRUN_RANDOM_186      unplaced-scaffold       na      na      KI270388.1      =       NT_187478.1     Primary Assembly        1216    chrUn_KI270388v1
        HSCHR19KIR_FH06_BA1_HAP_CTG3_1  alt-scaffold    19      Chromosome      KI270923.1      =       NT_187677.1     ALT_REF_LOCI_29 189352  chr19_KI270923v1_alt

The assembly report is found here: 
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/ or using [ftp provided by NCBI](ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.26.assembly.txt). 

Watch out for new versions. 

## Links 

* [Broad Bundle for hg38 at Google Cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)  
* [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml)
* [Link to pipeline resources by NCBI (GRCh38)](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)

Contains precompiled bwa, hisat2 and bowtie indexes with/without alt chromosomes and decoy chromosome, as well as RefSeq annotations.    

* [Archived data to match specific versions from ensembl](http://dec2017.archive.ensembl.org/Homo_sapiens/Info/Index)

Choose which patch you're at in the drop-down menu at the top. 


* [dbSNP with UCSC naming (in spite of source...)](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)
  All and common variants. [Human SNPs](ftp://ftp.ncbi.nih.gov/snp/organisms)  

* [dbSNP with NCBI naming](ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF) 
  

## hg38 UCSC genome 
Contains decoy sequences, HLA alleles,  
        
Example of fasta-headers from genome sequence: 
        
        >chr21  AC:CM000683.2  gi:568336003  LN:46709983  rl:Chromosome  M5:974dc7aec0b755b19f031418fdedf293  AS:GRCh38  hm:multiple
        >chr22_KI270733v1_random  AC:KI270733.1  gi:568335380  LN:179772  rg:chr22  rl:unlocalized  M5:f1fa05d48bb0c1f87237a28b66f0be0b  AS:GRCh38
        >chrUn_KI270753v1  AC:KI270753.1  gi:568335260  LN:62944  rl:unplaced  M5:25075fb2a1ecada67c0eb2f1fe0c7ec9  AS:GRCh38
        >chr4_KI270925v1_alt  AC:KI270925.1  gi:568335977  LN:555799  rg:chr4:189361858-189907070  rl:alt-scaffold  M5:08b3f5ff9beac37e3bbc50767cf3b43b  AS:GRCh38
        >HLA-DRB1*15:03:01:01   HLA00870 11567 bp
        >chrUn_JTFH01001740v1_decoy  AC:JTFH01001740.1  gi:725022526  LN:9344  rl:unplaced  M5:7c7fa5c413301e5c71020ee588c0c1cf  AS:hs38d1
        
        
       
## Discussions for tools 

* [BaseRecalibration](https://gatkforums.broadinstitute.org/gatk/discussion/6800/known-sites-for-indel-realignment-and-bqsr-in-hg38-bundle) 
 



## Resources for tools / steps
[There's the Broad Bundle.](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/) 
  
* Variant Callings (Mutect2, HaplotypeCaller (GATK) / Haplotyper (Sentieon)):
    * dbSNP, COSMIC, [gnomAD (giant file of 235 gb)](https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.s), [the smaller version of gnomAD](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/)
        * https://gnomad.broadinstitute.org/downloads 
    
  
     