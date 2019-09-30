# Release notes for pipeline versions 

## PST02 
_Dated 29-09-2019_ 
* Alignment  method, STAR, was updated to the newest version (2.7.2b) which facilitated getting counts for each gene to use in gene expression analysis (this method matches the output of htseq-count). 
* The new method also made it possible to remove samtools steps for sorting and adding read groups as the newest version of STAR also can take care of this without problems with memory. This almost removes a step in the pipeline, as well as writing a new sorted bam-file and a bam-file with readgroups. (We still need samtools to index the bam-files produced by STAR but this takes seconds.) 


## PST01
_Dated 20-08-2019_
* First pipeline version. 
