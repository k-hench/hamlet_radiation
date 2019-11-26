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

process het_inds {
	publishDir "../../2_analysis/het/raw", mode: 'copy'
	label "L_190g4h_inds"

	input:
	val( ind ) from inds_ch.splitCsv().map{it[0]}

	output:
	file( "*.het.gz" ) into inds_out

	script:
	"""
	vcftools --gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/filterd_bi-allelic.allBP.vcf.gz \
		--indv ${ind} \
		--counts2 \
		--stdout | \
		sed  's/{COUNT}/REF\\tALT/' | \
		cut -f 1,2,5,6 | \
		awk -v OFS='\\t' -v ind=${ind} '{if(NR==1){print \$1,\$2,"HET","IND"}else{x = \$3%2; ; print \$1,\$2,x/2,ind}}' | \
		gzip > ${ind}.het.gz
	"""
}

process win_inds {
	publishDir "../../2_analysis/het/50kb", mode: 'copy'
	label "L_20g2h_inds"
	module "R3.5.2"

	input:
	file( het ) from inds_out

	output:
	file( "*_win_het.gz" ) into win_out

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/het_by_ind.sh ${het} 50000
	"""
}
