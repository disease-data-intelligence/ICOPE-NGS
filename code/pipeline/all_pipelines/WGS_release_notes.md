# Release notes for pipeline versions 

_All minor updates are listed but do not affect the results only logistics and runtime._ 

## PSG01
_Dated 20-08-2019_
* First pipeline version.

### Minor updates

* _14-10-2019_ New naming of output vcf-files (independent of results) The main vcf-output created with Haplotyper is no longer called _-hc_, but just the sample name with vcf-suffix. 
* _24-11-2019_ New naming: Added prefix with sample identifier to all output files including run.log and quality_reports such that everything is named <sample>.quality_reports, <sample>.run.log etc. This also includes copy of the script that has been run (sentieon_WGS_PSG01.sh -> <sample>.sentieon_WGS_PSG01.sh)  
* _24-11-2019_ Added intermediate check points that will assess different number of quality parameters. If they are not passed, the pipeline will exit. 
The checks are of four kinds of checks and are relatively simple: 
    1. Checks number of paired reads in alignment and number of aligned paired reads.
    2. The same as second, but using the output from the run.log in different steps - for instance, if remove duplicates cuts away too many reads, the pipeline will discontinue. 
    3. Finally, the pipeline will asses the number of variants by counting the number of variant lines. 