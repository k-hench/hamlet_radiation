#!/usr/bin/env nextflow

// ----------------------- DISCLAIMER ----------------------
// this pipeline was not actually run using nextflow,
// but managed manually
// ---------------------------------------------------------

// git 17.1
// load genotypes
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ genotypes_ch }

// git 17.2
// drop serranus
process drop_serranus {

	input:
	set vcfId, file( vcf ) from genotypes_ch

	output:
	file( "hyp_mac2.vcf.gz" ) into no_serranus_ch

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep "tor\\|tab\\" > serr.pop

	vcftools --gzvcf ${vcf[0]} \
		--remove serr.pop \
		--recode \
		--stdout | bgzip > hyp_mac2.vcf.gz
	"""
}

// git 17.4
// LD filtering
process ld_filter {

	input:
	file( vcf ) from no_serranus_ch

	output:
	file( "hyp_mac2_ld05.vcf.gz" ) into ld_filtered_ch

	script:
	"""
	bcftools +prune \
		-l 0.5 \
		-w 50kb \
		${vcf} \
		-Oz \
		-o hyp_mac2_ld05.vcf.gz
	"""
}

// git 17.5
// load trios config
Channel
	.fromPath("../../ressources/hyp_sets.txt")
	.set{ hyp_sets_ch }

// git 17.6
// run D trios
process run_dtrios {
	publishDir "../../2_analysis/dstats/", mode: 'copy'

	input:
	set file( vcf ), file( sets ) from ld_filtered_ch.combine( hyp_sets_ch )

	output:
	set file( "hyp_ld05_dtrios_BBAA.txt" ), file( "hyp_ld05_dtrios_Dmin.txt" ) into dtrios_results_ch

	script:
	"""
	zcat ${vcf} > hyp_mac2_ld05.vcf

	Dsuite Dtrios \
		-c \
		-o hyp_ld05_dtrios \
		hyp_mac2_ld05.vcf \
		${sets}
	"""
}

// git 17.7
// load predefined species order
Channel
	.fromPath("../../ressources/species_order_alpha.txt")
	.set{ spec_order_ch }

// git 17.8
// Extract significant Dmin, BBAA trios
process run_correction {
	publishDir "../../2_analysis/dstats/", mode: 'copy'

	input:
	set file( vcf ), file( sets ), file( spec_order) from dtrios_results_ch.combine( spec_order_ch )

	output:
	set file( "BBAA_sign_ld05.csv" ), file( "Dmin_sign_ld05.csv" ) into dtrios_signif_ch

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/dstats.R
	"""
}
