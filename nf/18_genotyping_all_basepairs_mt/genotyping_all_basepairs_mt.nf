#!/usr/bin/env nextflow
// git 18.1
// open genotype likelyhoods
Channel
	.fromFilePairs("../../1_genotyping/1_gvcfs/cohort.g.vcf.{gz,gz.tbi}")
	.set{ vcf_cohort }

Channel
	.from(["LG_M", "unplaced"])
	.set{ lg_mode }

// git 18.2
// actual genotyping step (including invariant sites)
process joint_genotype_snps {
	label "L_88g48h_LGs_genotype"
	publishDir "../../1_genotyping/2_raw_vcfs/", mode: 'copy'

	input:
	set vcfId, file( vcf ), val( mode ) from vcf_cohort.combine( lg_mode )

	output:
	set file( "all_sites.${mode}.vcf.gz" ), file( "all_sites.${mode}.vcf.gz.tbi" ), val( mode ) into ( all_bp_non_lg_1, all_bp_non_lg_2 )

	script:
	if( mode == 'unplaced' )
	"""
	gatk --java-options "-Xmx85g" \
		GenotypeGVCFs \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-XL=LG01 \
		-XL=LG02 \
		-XL=LG03 \
		-XL=LG04 \
		-XL=LG05 \
		-XL=LG06 \
		-XL=LG07 \
		-XL=LG08 \
		-XL=LG09 \
		-XL=LG10 \
		-XL=LG11 \
		-XL=LG12 \
		-XL=LG13 \
		-XL=LG14 \
		-XL=LG15 \
		-XL=LG16 \
		-XL=LG17 \
		-XL=LG18 \
		-XL=LG19 \
		-XL=LG20 \
		-XL=LG21 \
		-XL=LG22 \
		-XL=LG23 \
		-XL=LG24 \
		-XL=LG_M \
		-V=${vcf[0]} \
		-O=intermediate.vcf.gz \
		--include-non-variant-sites=true \
		--allow-old-rms-mapping-quality-annotation-data

	gatk --java-options "-Xmx85G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		--select-type-to-exclude=INDEL \
		-O=all_sites.${mode}.vcf.gz

	rm intermediate.*
	"""
	else if( mode == 'LG_M' )
	"""
	gatk --java-options "-Xmx85g" \
		GenotypeGVCFs \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-L=${mode} \
		-V=${vcf[0]} \
		-O=intermediate.vcf.gz \
		--include-non-variant-sites=true \
		--allow-old-rms-mapping-quality-annotation-data

	gatk --java-options "-Xmx85G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		--select-type-to-exclude=INDEL \
		-O=all_sites.${mode}.vcf.gz

	rm intermediate.*
	"""
}
