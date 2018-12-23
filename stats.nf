#!/usr/bin/env nextflow

Channel
	.fromFilePairs('\$BASE_DIR/1_genotyping/4_phased/phased_mac2.vcf.{gz,.gz.tbi}')
	.into{ vcf_phylo; vcf_locations }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

process subset_vcf_by_location {
	   label "L_20g2h_subset_vcf"

	   input:
		 val( loc ) from locations_ch
		 set vcfId, file( vcf ) from vcf_locations

	   output:
	   set val( loc ) file( "${loc}.vcf.gz" ) into ( vcf_loc_pca, vcf_loc_fst, vcf_loc_admix )

	   script:
	   """
		 awk v loc="${loc}" '{if( ( !($3=="tor" || $3=="tab" ) && NR>1 ) && $4==loc ){print}}' \$BASE_DIR/metadata/sample_info.txt | \
		 	cut -f 1 > ${loc}.pop

	   vcftools --gzvcf ${vcf[0]} \
	    --keep ${loc}.pop \
	    --mac 1 \
	    --recode \
	    --stdout | bgzip > ${loc}.vcf.gz
	   """
	 }

/*
"""
Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${loc}.vcf.gz \$BASE_DIR/metadata/sample_info.txt 8
"""
*/