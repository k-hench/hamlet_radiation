#!/usr/bin/env nextflow
/* create channel of linkage groups */
Channel
	.from( 2..12 )
	.into{ admx_ch; admx_loc_ch }

Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations; vcf_all_samples_pca; vcf_admx; vcf_geno }

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
	set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_pca, vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3, vcf_loc_admix )

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep ${loc} | \
		grep -v tor | \
		grep -v tab > ${loc}.pop
	vcftools --gzvcf ${vcf[0]} \
		--keep ${loc}.pop \
		--mac 3 \
		--recode \
		--stdout | gzip > ${loc}.vcf.gz
	"""
}

/* 1) PCA section ============== */
/* 1a) PCA (local) -------------- */
process pca_location {
	label "L_20g15h_pca_location"
	publishDir "../../figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "../../2_analysis/pca", mode: 'copy' , pattern: "*.gz"

	input:
	set val( loc ), file( vcf ), file( pop ) from vcf_loc_pca

	output:
	set file( "${loc}.prime_pca.pdf" ), file( "${loc}.pca.pdf" ), file( "${loc}.exp_var.txt.gz" ), file( "${loc}.scores.txt.gz" ) into pca_loc_out

	script:
	"""
	awk '{print \$1"\\t"\$1}' ${loc}.pop | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > ${loc}.pop.txt
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R ${loc}.pop.txt 6
	"""
}

/* 1b) PCA (global) -------------- */
process pca_all {
	label "L_20g15h_pca_all"
	publishDir "../../figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "../../2_analysis/pca", mode: 'copy' , pattern: "*.txt.gz"
	publishDir "../../1_genotyping/4_phased/", mode: 'copy' , pattern: "*.vcf.gz"

	input:
	set vcfId, file( vcf ) from vcf_all_samples_pca

	output:
	set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out
	file( "hamlets_only.vcf.gz*" ) into vcf_hamlets_only
	set file( "hamlets_only.vcf.gz*" ), file( "hamlets_only.pop.txt" ) into vcf_multi_fst

	script:
	"""
	# complete PCA, all samples ------------
	vcfsamplenames ${vcf[0]} | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > all.pop.txt
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R all.pop.txt 6
	# PCA without outgroups ---------------
	vcfsamplenames ${vcf[0]} | \
		grep -v "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" > outgroup.pop
	vcfsamplenames ${vcf[0]} | \
		grep "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > hamlets_only.pop.txt
	vcftools \
		--gzvcf ${vcf[0]} \
		--remove outgroup.pop \
		--recode \
		--stdout | gzip > hamlets_only.vcf.gz
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.vcf.gz \$BASE_DIR/R/project_config.R hamlets_only.pop.txt 6
	# PCA without indtgo or gumigutta ---------------
	grep -v "ind\\|gum" hamlets_only.pop.txt > core_hamlets.pop.txt
	cut -f 1 core_hamlets.pop.txt > core_hamlets.pop
	vcftools \
		--gzvcf ${vcf[0]} \
		--keep core_hamlets.pop \
		--recode \
		--stdout | gzip > core_hamlets.vcf.gz
	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R core_hamlets.vcf.gz \$BASE_DIR/R/project_config.R core_hamlets.pop.txt 6
	"""
}
