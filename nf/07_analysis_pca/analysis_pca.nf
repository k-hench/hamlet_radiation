#!/usr/bin/env nextflow
// git 7.1
// prepare subset modes (whole genome vs non-diverged regions)
Channel
	.from( "whg", "subset_non_diverged")
	.set{ subset_type_ch }

// git 7.2
// load table with differentiation outlier regions
Channel
	.fromPath( "../../2_analysis/summaries/fst_outliers_998.tsv" )
	.set{ outlier_tab }

// git 7.3
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.combine( outlier_tab )
	.combine( subset_type_ch )
	.set{ vcf_ch }

// git 7.4
// depending on subset mode, subset vcf
process subset_vcf_divergence_based {
	label "L_20g2h_subset_divergence"

	input:
	set  vcfId, file( vcf ), file( outlier_tab ), val( subset_type ) from vcf_ch

	output:
	set file( "${subset_type}.vcf.gz" ), file( "${subset_type}.vcf.gz.tbi" ), val( subset_type ) into ( vcf_locations, vcf_all_samples_pca )

	script:
	"""
	if [ "${subset_type}" == "subset_non_diverged" ];then
		awk -v OFS="\\t" '{print \$2,\$3,\$4}' ${outlier_tab} > diverged_regions.bed 
		SUBSET="--exclude-bed diverged_regions.bed"
	else
		SUBSET=""
	fi

	vcftools --gzvcf ${vcf[0]} \
		\$SUBSET \
		--recode \
		--stdout | bgzip > ${subset_type}.vcf.gz
	
	tabix ${subset_type}.vcf.gz
	"""
}

// git 7.5
// prepare location channel for separate pcas
Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

// git 7.6
// attach genotypes to location channel
locations_ch
	.combine( vcf_locations )
	.set{ vcf_location_combo }

// git 7.7
// subset vcf by location
process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), file( vcf ), file( vcfidx ), val( subset_type ) from vcf_location_combo

	output:
	set val( loc ), file( "*.vcf.gz" ), file( "*.pop" ), val( subset_type ) into ( vcf_loc_pca )

	script:
	"""
	vcfsamplenames ${vcf} | \
		grep ${loc} | \
		grep -v tor | \
		grep -v tab > ${loc}.${subset_type}.pop
	vcftools --gzvcf ${vcf} \
		--keep ${loc}.${subset_type}.pop \
		--mac 3 \
		--recode \
		--stdout | gzip > ${loc}.${subset_type}.vcf.gz
	"""
}

// PCA section
// -----------
// git 7.8
// run pca by location
process pca_location {
	label "L_20g15h_pca_location"
	publishDir "../../figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "../../2_analysis/pca", mode: 'copy' , pattern: "*.gz"

	input:
	set val( loc ), file( vcf ), file( pop ), val( subset_type ) from vcf_loc_pca

	output:
	set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_loc_out

	script:
	"""
	awk '{print \$1"\\t"\$1}' ${loc}.${subset_type}.pop | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > ${loc}.${subset_type}.pop.txt
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf} ${loc}.${subset_type}.pop.txt 6
	"""
}

// git 7.9
// run pca for global data set
process pca_all {
	label "L_20g15h_pca_all"
	publishDir "../../figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "../../2_analysis/pca", mode: 'copy' , pattern: "*.txt.gz"
	publishDir "../../1_genotyping/4_phased/", mode: 'copy' , pattern: "*.vcf.gz"

	input:
	set file( vcf ), file( vcfidx ), val( subset_type ) from vcf_all_samples_pca

	output:
	set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out

	script:
	"""
	# complete PCA, all samples ------------
	vcfsamplenames ${vcf} | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > all.${subset_type}.pop.txt
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf} all.${subset_type}.pop.txt 6

	# PCA without outgroups ---------------
	vcfsamplenames ${vcf} | \
		grep -v "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni\\|flo" > outgroup.${subset_type}.pop
	vcfsamplenames ${vcf} | \
		grep "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni\\|flo" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > hamlets_only.${subset_type}.pop.txt
	vcftools \
		--gzvcf ${vcf} \
		--remove outgroup.${subset_type}.pop \
		--recode \
		--stdout | gzip > hamlets_only.${subset_type}.vcf.gz
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.${subset_type}.vcf.gz hamlets_only.${subset_type}.pop.txt 6
	"""
}
