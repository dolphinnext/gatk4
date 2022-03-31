Gatk4 pipeline is adopted from;
https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

This pipeline is intended for calling variants in samples that are clonal – i.e. a single individual.  
The frequencies of variants in these samples are expected to be 1 (for haploids or homozygous diploids) 
or 0.5 (for heterozygous diploids).  

Base Quality Score Recalibration (BQSR) is an important step for accurate variant detection that aims 
to minimize the effect of technical variation on base quality scores (measured as Phred scores). 