# Release notes for pipeline versions 

_All minor updates are listed but do not affect the results only logistics and runtime._  
## PST02 
_Dated 29-09-2019_ 
* Alignment  method, STAR, was updated to the newest version (2.7.2b) which facilitated getting counts for each gene to use in gene expression analysis (this method matches the output of htseq-count). 
* The new method also made it possible to remove samtools steps for sorting and adding read groups as the newest version of STAR also can take care of this without problems with memory. This almost removes a step in the pipeline, as well as writing a new sorted bam-file and a bam-file with readgroups. (We still need samtools to index the bam-files produced by STAR but this takes seconds.) 
### Minor updates

* _14-10-2019_ New naming of output vcf-files. The main vcf-output created with Haplotyper is no longer called _-hc_, but just the sample name with vcf-suffix. The vcf-file made with the specific tumour variantcaller is still kept but has the suffix _-TNScope_.     

## PST01
_Dated 20-08-2019_
* First pipeline version. 
