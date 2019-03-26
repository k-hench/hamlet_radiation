#!/usr/bin/env nextflow
Channel
	.fromFilePairs("../../1_genotyping/1_gvcfs/cohort.g.vcf.{gz,gz.tbi}")
	.set{ vcf_cohort }

Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ ch_LG_ids }

ch_LG_ids.combine( vcf_cohort ).set{ vcf_lg_combo }

/* actual genotyping step (varinat sites only) */
process joint_genotype_snps {
	label "L_O88g90h_LGs_genotype"

	input:
	set val( lg ), vcfId, file( vcf ) from vcf_lg_combo

	output:
	set val( 'all' ), val( lg ), file( "all_site*.vcf.gz" ), file( "all_site*.vcf.gz.tbi" ) into all_bp_by_location

	script:
	"""
	gatk --java-options "-Xmx85g" \
		GenotypeGVCFs \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-L=LG${lg} \
		-V=${vcf[0]} \
		-O=intermediate.vcf.gz \
		--include-non-variant-sites=true

	gatk --java-options "-Xmx85G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		--select-type-to-exclude=INDEL \
		-O=all_sites.LG${lg}.vcf.gz

	rm intermediate.*
	"""
}

process merge_genotypes {
	label 'L_78g5h_merge_genotypes'
	echo true

	input:
	set val( dummy ),  val( lg ), file( vcf ), file( tbi ) from all_bp_by_location.groupTuple()

	output:
	file( "all_sites.vcf.gz" ) into all_bp_merged

	script:
	"""
	INPUT=\$(ls -1 *vcf.gz | sed 's/^/ -I /g' | cat \$( echo ))

	gatk --java-options "-Xmx85g" \
		GatherVcfs \
		\$INPUT \
		-O=all_sites.vcf.gz
	"""
}

process filterSNPs {
	label 'L_105g30h_filter_genotypes'
	publishDir "../../1_genotyping/3_gatk_filtered/", mode: 'copy'
	module "openssl1.0.2"

	input:
	file( vcf ) from all_bp_merged

	output:
	set file( "filterd_bi-allelic.vcf.gz" ), file( "filterd_bi-allelic.vcf.gz.tbi" ) into filtered_snps

	script:
	"""
	tabix -p vcf ${vcf}

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
		--filter-name "filter_ReadPosRankSum"

	gatk --java-options "-Xmx75G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		-O=intermediate.filterd.vcf.gz \
		--exclude-filtered

	vcftools \
		--gzvcf intermediate.filterd.vcf.gz \
		--max-missing-count 17 \
		--max-alleles 2 \
		--stdout  \
		--recode | \
		bgzip > filterd_bi-allelic.allBP.vcf.gz

	tabix -p vcf filterd_bi-allelic.allBP.vcf.gz

	rm intermediate.*
	"""
}
