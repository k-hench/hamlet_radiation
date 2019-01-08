#!/usr/bin/env nextflow

Channel
	.fromFilePairs('\$BASE_DIR/1_genotyping/4_phased/phased_mac2.vcf.{gz,.gz.tbi}')
	.into{ vcf_phylo; vcf_locations; vcf_all_samples_pca }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

vcf_location_combo = locations_ch.combine( vcf_locations )

process subset_vcf_by_location {
	   label "L_20g2h_subset_vcf"

	   input:
		 set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	   output:
	   set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_pca, vcf_loc_fst, vcf_loc_admix )

	   script:
	   """
			vcfsamplenames phased_mac2.vcf.gz | \
				grep ${loc} | \
				grep -v tor | \
				grep -v tab > ${loc}.pop

	   vcftools --gzvcf ${vcf[0]} \
	    --keep ${loc}.pop \
	    --mac 1 \
	    --recode \
	    --stdout | bgzip > ${loc}.vcf.gz
	   """
	 }
/*
	process pca_location {
			label "L_20g15h_pca_location"
			publishDir "figures/pca", mode: 'move'

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

process pca_all {
		label "L_20g15h_pca_all"
		publishDir "figures/pca", mode: 'move'

		input:
		file( vcf ) from vcf_all_samples_pca

		output:
		set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out

		script:
		"""
		vcfsamplenames ${vcf[0]} | \
			'{print \$1"\\t"\$1}' | \
			sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > all.pop.txt

		Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R all.pop.txt 6
		"""
}
*/
