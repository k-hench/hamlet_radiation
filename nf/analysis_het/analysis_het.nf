#!/usr/bin/env nextflow
// nextflow run het.nf -c nextflow.config
// git 7.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_by_ind; }

process split_inds {
	label "L_loc_split_vcf"

	input:
	set vcfId, file( vcf ) from vcf_by_ind

	output:
	file( "inds.txt" ) into inds_ch

	script:
	"""
	vcfsamplenames ${vcf[0]} > inds.txt
	"""
}

process get_inds {
	publishDir "../../2_analysis/het", mode: 'copy'
	label "L_20g2h_inds"

	input:
	val( ind ) from inds_ch.splitCsv().map{it[0]}

	output:
	file( "*.het.gz" ) into inds_out

	script:
	"""
	vcftools --gzvcf \$BASE_DIR/1_genotyping/4_phased/phased_mac2.vcf.gz \
		--indv ${ind} \
		--counts2 \
		--stdout | \
		sed  's/{COUNT}/REF\\tALT/' | \
		cut -f 1,2,5,6 | \
		awk -v OFS='\\t' -v ind=${ind} '{if(NR==1){print \$1,\$2,"HET","IND"}else{x = \$3%2; ; print \$1,\$2,x/2,ind}}' | \
		gzip > ${ind}.het.gz
	"""
}