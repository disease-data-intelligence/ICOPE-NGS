# Description 
Pipeline for paired somatic variant calling. Usually used with the submission script `somatic_setup.py`. 
It needs a tumor and germline bam-file and a destination path (for a folder), which is a combination of the two. 
For instance: 

```
# example of submission 
germline=11021/11021_154854.G.XXX/11021_154854.G.XXX.SX/11021_154854.G.XXX.SX.PSG01/11021_154854.G.XXX.SX.PSG01.bam
tumor=11021/11021_154854.T.XXX/11021_154854.T.XXX.SX/11021_154854.T.XXX.SX.PST02/11021_154854.T.XXX.SX.PST02.bam
destination=11021/11021_154854.T.XXX/11021_154854.T.XXX.SX/11021-154854.T.XXX.SX.PST02-154854.G.XXX.SX.PSG01-PSP01
sentieon_paired.sh $germline $tumor $destination
```
Meaning that the destination folder is `<MRD>-<full_tumor_name>-<full_germline_name>-PSP01` and is placed in the tumor sequencing folder. 

# Release notes for pipeline versions 

_All minor updates are listed but do not affect the results only logistics and runtime._ 

## PSP01
_Dated 30-09-2019_
* First paired pipeline version for finding somatic variants.

### Minor updates

* _14-10-2019_ New naming of output vcf-files. The main somatic variant output is no longer called  _-TNScope_, but just the sample name with vcf-suffix.  
