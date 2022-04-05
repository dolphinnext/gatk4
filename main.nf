$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.bwaref){params.bwaref = ""} 
if (!params.firstRound){params.firstRound = ""} 
if (!params.secondRound){params.secondRound = ""} 
if (!params.snpeff_db){params.snpeff_db = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_1_reads_g_0}
 } else {  
	g_1_reads_g_0 = Channel.empty()
 }

Channel.value(params.mate).set{g_2_mate_g_0}
Channel.value(params.bwaref).into{g_5_ref_flat_g_0;g_5_ref_flat_g_4;g_5_ref_flat_g_6;g_5_ref_flat_g_20;g_5_ref_flat_g_8;g_5_ref_flat_g_22;g_5_ref_flat_g_9;g_5_ref_flat_g_24;g_5_ref_flat_g_14;g_5_ref_flat_g_16}
Channel.value(params.firstRound).set{g_7_roundBatchAlgorithm_g_6}
Channel.value(params.secondRound).set{g_21_roundBatchAlgorithm_g_20}
Channel.value(params.snpeff_db).set{g_28_ref_flat_g_26}


process align {

input:
 set val(name),file(reads) from g_1_reads_g_0
 val mate from g_2_mate_g_0
 val ref from g_5_ref_flat_g_0

output:
 set val(name), file("${name}_aligned_reads.bam")  into g_0_bam_file0_g_35

script:

readGroup = "@RG\\tID:${name}\\tLB:${name}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${name}"

"""
    bwa mem \
	-K 100000000 \
	-v 3 \
	-t 1 \
	-Y \
	-R '${readGroup}' \
	${ref} \
	${reads} \
	> ${name}_aligned_reads.sam
	
	samtools view -bS ${name}_aligned_reads.sam > ${name}_aligned_reads.bam 
"""
}


process markDuplicatesSpark {

input:
 set val(name), file(reads) from g_0_bam_file0_g_35

output:
 set val(name), file("${name}_sorted_dedup.bam")  into g_35_mapped_reads0_g_4, g_35_mapped_reads0_g_6
 file "${name}_dedup_metrics.txt"  into g_35_txtFile11

script:
"""
. /opt/conda/etc/profile.d/conda.sh
conda activate base
mkdir -p tmp/${name}
gatk --java-options '-Djava.io.tmpdir=tmp/${name}' \
 MarkDuplicatesSpark \
-I ${reads} \
-M ${name}_dedup_metrics.txt \
-O ${name}_sorted_dedup.bam 
rm -r tmp/${name}
"""
}


process HaplotypeCaller {

input:
 set val(name), file(input_bam) from g_35_mapped_reads0_g_6
 val ref from g_5_ref_flat_g_6
 val round from g_7_roundBatchAlgorithm_g_6

output:
 set val(name), val(round), file("${input_bam}"), file("${name}_raw_variants_${round}.vcf")  into g_6_VCFset0_g_8

script:
"""
 gatk HaplotypeCaller \
	-R $ref \
	-I $input_bam \
	-O ${name}_raw_variants_${round}.vcf
"""
}


process selectVariants {

input:
 set val(name), val(round),  file(inputbam), file(variants) from g_6_VCFset0_g_8
 val ref from g_5_ref_flat_g_8

output:
 set val(name), val(round),  file("${inputbam}"), file("${name}_raw_*_${round}.vcf*")  into g_8_VCFset0_g_9

script:
"""
gatk SelectVariants \
	-R ${ref} \
	-V ${variants} \
	-select-type SNP \
	-O ${name}_raw_snps_${round}.vcf
gatk SelectVariants \
	-R ${ref} \
	-V ${variants} \
	-select-type INDEL \
	-O ${name}_raw_indels_${round}.vcf
"""
}


process VariantFiltration {

input:
 set val(name), val(round), file(inputbam), file(variants) from g_8_VCFset0_g_9
 val ref from g_5_ref_flat_g_9

output:
 set val(name), val(round), file("${inputbam}"), file("${name}_filtered_*_${round}.vcf*")  into g_9_VCFset0_g_14

script:
	
"""
#FOR SNPS
gatk VariantFiltration \
	-R ${ref} \
	-V ${name}_raw_snps_${round}.vcf \
	-O ${name}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#FOR INDELS
gatk VariantFiltration \
        -R ${ref} \
        -V ${name}_raw_indels_${round}.vcf \
        -O ${name}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"

"""
}


process BaseRecalibrator {

input:
 set val(name), val(round), file(input_bam), file(filtered_variants) from g_9_VCFset0_g_14
 val ref from g_5_ref_flat_g_14

output:
 set val(name), file("${input_bam}"), file("${name}_bqsr_*.vcf*"), file("${name}_recal_data.table")  into g_14_VCFset0_g_16

script:
"""
gatk SelectVariants \
	--exclude-filtered \
	-V ${name}_filtered_snps_${round}.vcf \
	-O ${name}_bqsr_snps.vcf

gatk SelectVariants \
	--exclude-filtered \
	-V ${name}_filtered_indels_${round}.vcf \
	-O ${name}_bqsr_indels.vcf
	
gatk BaseRecalibrator \
	-R ${ref} \
	-I ${input_bam} \
	--known-sites ${name}_bqsr_snps.vcf \
    --known-sites ${name}_bqsr_indels.vcf \
	-O ${name}_recal_data.table

"""
}


process applyBSQRS {

input:
 set val(name), file(input_bam), file(variants), file(recal_table) from g_14_VCFset0_g_16
 val ref from g_5_ref_flat_g_16

output:
 set val(name),  file("${name}_recal_data.table"),  file("${name}_post_recal_data.table")  into g_16_outputFileTab0_g_17
 set val(name), file("${name}_recal.bam")  into g_16_mapped_reads1_g_20

script:
"""
    . /opt/conda/etc/profile.d/conda.sh
    conda activate base
	gatk ApplyBQSR \
        -R $ref \
        -I $input_bam \
        -bqsr ${name}_recal_data.table \
        -O ${name}_recal.bam
    gatk BaseRecalibrator \
	    -R ${ref} \
		-I ${name}_recal.bam \
	    --known-sites ${name}_bqsr_snps.vcf\
		--known-sites ${name}_bqsr_indels.vcf \
		-O ${name}_post_recal_data.table
"""
}


process AnalyzeCovariates {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_recalibration_plots.pdf$/) "AnalyzeCovarage/$filename"}
input:
 set val(name),file(recal_table), file(post_recal_table) from g_16_outputFileTab0_g_17

output:
 file "${name}_recalibration_plots.pdf"  into g_17_outputFilePdf00

script:

"""
gatk AnalyzeCovariates \
	-before $recal_table \
	-after $post_recal_table \
	-plots ${name}_recalibration_plots.pdf
"""
}


process getMetrics {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_insert_size_histogram.pdf$/) "report/$filename"}
input:
 set val(name),file(sorted_dedup_reads) from g_35_mapped_reads0_g_4
 val ref from g_5_ref_flat_g_4

output:
 file "${name}_alignment_metrics.txt"  into g_4_txtFile00
 file "${name}_insert_metrics.txt"  into g_4_txtFile11
 file "${name}_insert_size_histogram.pdf"  into g_4_outputFilePdf22
 file "${name}_depth_out.txt"  into g_4_txtFile33

errorStrategy 'retry'
maxRetries 1

when:
(params.run_Metrics && (params.run_Metrics == "yes")) || !params.run_Metrics

script:
"""
    picard \
       CollectAlignmentSummaryMetrics \
	   R=${ref} \
       I=${sorted_dedup_reads} \
	   O=${name}_alignment_metrics.txt
    picard \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	    OUTPUT=${name}_insert_metrics.txt \
        HISTOGRAM_FILE=${name}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${name}_depth_out.txt
"""
}


process CalibHaplotypeCaller {

input:
 set val(name), file(input_bam) from g_16_mapped_reads1_g_20
 val ref from g_5_ref_flat_g_20
 val round from g_21_roundBatchAlgorithm_g_20

output:
 set val(name), val(round), file("${input_bam}"), file("${name}_raw_variants_${round}.vcf")  into g_20_VCFset0_g_22

script:
"""
 gatk HaplotypeCaller \
	-R $ref \
	-I $input_bam \
	-O ${name}_raw_variants_${round}.vcf
"""
}


process selectPostSNPs {

input:
 set val(name), val(round),  file(inputbam), file(variants) from g_20_VCFset0_g_22
 val ref from g_5_ref_flat_g_22

output:
 set val(name), val(round),  file("${inputbam}"), file("${name}_raw_*_${round}.vcf*")  into g_22_VCFset0_g_24

script:
"""
gatk SelectVariants \
	-R ${ref} \
	-V ${variants} \
	-select-type SNP \
	-O ${name}_raw_snps_${round}.vcf
gatk SelectVariants \
	-R ${ref} \
	-V ${variants} \
	-select-type INDEL \
	-O ${name}_raw_indels_${round}.vcf
"""
}


process SNPPostFIltration {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /(${inputbam}|${name}_filtered_.*_${round}.vcf.*)$/) "postfilter_vcfout/$filename"}
input:
 set val(name), val(round), file(inputbam), file(variants) from g_22_VCFset0_g_24
 val ref from g_5_ref_flat_g_24

output:
 set val(name), val(round), file("${inputbam}"), file("${name}_filtered_*_${round}.vcf*")  into g_24_VCFset0_g_26

script:
	
"""
#FOR SNPS
gatk VariantFiltration \
	-R ${ref} \
	-V ${name}_raw_snps_${round}.vcf \
	-O ${name}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#FOR INDELS
gatk VariantFiltration \
        -R ${ref} \
        -V ${name}_raw_indels_${round}.vcf \
        -O ${name}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"

"""
}


process SnpEff {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtered_snps.ann.vcf$/) "snp_vcfout/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_snpEff_summary.html$/) "SnpEff_report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_snpEff_genes.txt$/) "SnpEff/$filename"}
input:
 set val(name), val(round), file(filtered_snps), file(filtered_index)  from g_24_VCFset0_g_26
 val snpeff_db from g_28_ref_flat_g_26

output:
 set val(name), file("${name}_filtered_snps.ann.vcf")  into g_26_VCFset00
 file "${name}_snpEff_summary.html"  into g_26_outputHTML11
 file "${name}_snpEff_genes.txt"  into g_26_txtFile22

script:
"""
# R64-1-1.86 Saccharomyces_cerevisiae
# GRCh38.p7.RefSeq Homo_sapiens
# GRCm38.75 Mus_musculus
. /opt/conda/etc/profile.d/conda.sh 
conda activate base
snpEff -v \
	-dataDir /tmp/data \
	${snpeff_db} \
	$filtered_snps > ${name}_filtered_snps.ann.vcf 2>out.log
	mv snpEff_summary.html ${name}_snpEff_summary.html
	sed -i -e '1d' snpEff_genes.txt 
	mv snpEff_genes.txt ${name}_snpEff_genes.txt
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
