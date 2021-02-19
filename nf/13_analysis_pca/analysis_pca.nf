#!/usr/bin/env nextflow
Channel
	.from( "whg", "subset_non_diverged")
	.set{ subset_type_ch }

Channel
	.fromPath( "../../2_analysis/summaries/fst_outliers_998.tsv" )
	.set{ outlier_tab }

Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.combine( outlier_tab )
	.combine( subset_type_ch )
	.set{ vcf_ch }

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
		SUBSET="--exclude-positions diverged_regions.bed"
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
//	.into{ vcf_locations; vcf_all_samples_pca}

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

locations_ch
	.combine( vcf_locations )
	.set{ vcf_location_combo }

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

// 1) PCA section ==============
// 1a) PCA (local) --------------
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

// 1b) PCA (global) -------------- 
process pca_all {
	label "L_20g15h_pca_all"
	publishDir "../../figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "../../2_analysis/pca", mode: 'copy' , pattern: "*.txt.gz"
	publishDir "../../1_genotyping/4_phased/", mode: 'copy' , pattern: "*.vcf.gz"

	input:
	set file( vcf ), file( vcfidx ), val( subset_type ) from vcf_all_samples_pca

	output:
	set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out
	file( "hamlets_only.${subset_type}.vcf.gz" ) into vcf_hamlets_only
	set file( "hamlets_only.${subset_type}.vcf.gz" ), file( "hamlets_only.${subset_type}.pop.txt" ) into vcf_multi_fst

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
