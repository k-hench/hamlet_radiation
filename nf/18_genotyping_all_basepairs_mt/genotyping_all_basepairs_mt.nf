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

// git 1.13
/* produce metrics table to determine filtering thresholds - ups forgot to extract SNPS first*/
process joint_genotype_metrics {
	label 'L_28g5h_genotype_metrics'
	publishDir "../../1_genotyping/2_raw_vcfs/", mode: 'copy'

	input:
	set file( vcf ), file( tbi ), val( mode )  from all_bp_non_lg_1

	output:
	file( "${vcf}.${mode}.table.txt" ) into raw_metrics

	script:
	"""
	gatk --java-options "-Xmx25G" \
		VariantsToTable \
		--variant=${vcf} \
		--output=${vcf}.${mode}.table.txt \
		-F=CHROM -F=POS -F=MQ \
		-F=QD -F=FS -F=MQRankSum -F=ReadPosRankSum \
		--show-filtered
	"""
}

/*
// git 18.3
// quality based filtering
process filterSNP_first {
	label 'L_88g48h_filter_gt1'

	input:
	set file( vcf ), file( tbi ), val( mode )  from all_bp_non_lg_2

	output:
	set file( "intermediate.filterd.${mode}.vcf.gz" ), file( "intermediate.filterd.${mode}.vcf.gz.tbi" ) into filtered_snps_first

	script:
	"""
	gatk --java-options "-Xmx75G" \
		VariantFiltration \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V ${vcf} \
		-O=intermediate.vcf.gz \
		--filter-expression "QD < 2.5" \
		--filter-name "filter_QD" \
		--filter-expression "FS > 25.0" \
		--filter-name "filter_FS" \
		--filter-expression "MQ < 52.0 || MQ > 65.0" \
		--filter-name "filter_MQ" \
		--filter-expression "MQRankSum < -0.2 || MQRankSum > 0.2" \
		--filter-name "filter_MQRankSum" \
		--filter-expression "ReadPosRankSum < -2.0 || ReadPosRankSum > 2.0 " \
		--filter-name "filter_ReadPosRankSum" \
		--QUIET true &> var_filt.log

	gatk --java-options "-Xmx75G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		-O=intermediate.filterd.${mode}.vcf.gz \
		--exclude-filtered \
		--QUIET true \
		--verbosity ERROR  &> var_select.log
	"""
}

// git 18.4
// missingness based filtering
// the resulting vcf file represents
// the 'all BP' data set
process filterSNP_second {
	label 'L_88g48h_filter_gt2'
	publishDir "../../1_genotyping/3_gatk_filtered/", mode: 'copy'

	input:
	set file( vcf ), file( tbi ) from filtered_snps_first

	output:
	file( "filterd.allBP.nonLG.vcf.gz" ) into filtered_snps

	script:
	"""
	vcftools \
		--gzvcf ${vcf} \
		--max-missing-count 17 \
		--stdout  \
		--recode | \
		gzip > filterd.allBP.nonLG.vcf.gz
	"""
}
*/
// alias nf_run_allbp_mt="nextflow run genotyping_all_basepairs_mt.nf -with-dag docs/genotyping_allBP_mt.dot -with-report ../../docs/genotyping_allBP_mt.html -with-timeline ../../docs/genotyping_allBP_mt_timeline.html -c ../../nextflow.config -resume"
// alias nf_run_allbp_mt="nextflow run genotyping_all_basepairs_mt.nf -c ../../nextflow.config -resume"
