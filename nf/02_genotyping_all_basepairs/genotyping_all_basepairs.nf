#!/usr/bin/env nextflow
// git 2.1
// open genotype likelyhoods
Channel
	.fromFilePairs("../../1_genotyping/1_gvcfs/cohort.g.vcf.{gz,gz.tbi}")
	.set{ vcf_cohort }

// git 2.2
// initialize LG channel
Channel
	.from(["_M", "unplaced"])
	.set{ lg_mode }

Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.concat( lg_mode )
	.set{ ch_LG_ids }

// git 2.3
// combine genotypes and LGs
ch_LG_ids.combine( vcf_cohort ).set{ vcf_lg_combo }

// git 2.4
// actual genotyping step (including invariant sites)
process joint_genotype_snps {
	label "L_O88g90h_LGs_genotype"

	input:
	set val( lg ), vcfId, file( vcf ) from vcf_lg_combo

	output:
	set val( 'all' ), val( lg ), file( "all_site*.vcf.gz" ), file( "all_site*.vcf.gz.tbi" ) into all_bp_by_location

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
		-O=all_sites.${lg}.vcf.gz

	rm intermediate.*
	"""
	else
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

// git 2.5
// merge all LGs
process merge_genotypes {
	label 'L_78g5h_merge_genotypes'
	echo true

	input:
	set val( dummy ), val( lg ), file( vcf ), file( tbi ) from all_bp_by_location.filter({ it[1] =~/^[0-2].*/ }).groupTuple()

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

// git 2.6
// quality based filtering
process filterSNP_first {
	label 'L_105g30h_filter_gt1'

	input:
	file( vcf ) from all_bp_merged

	output:
	set file( "intermediate.filterd.vcf.gz" ), file( "intermediate.filterd.vcf.gz.tbi" ) into filtered_snps_first

	script:
	"""
	module load openssl1.0.2

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
		--filter-name "filter_ReadPosRankSum" \
		--QUIET true &> var_filt.log

	gatk --java-options "-Xmx75G" \
		SelectVariants \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V=intermediate.vcf.gz \
		-O=intermediate.filterd.vcf.gz \
		--exclude-filtered \
		--QUIET true \
		--verbosity ERROR  &> var_select.log
	"""
}

// git 2.7
// missingness based filtering
// the resulting vcf file represents
// the 'all BP' data set
process filterSNP_second {
	label 'L_105g30h_filter_gt2'
	publishDir "../../1_genotyping/3_gatk_filtered/", mode: 'copy'

	input:
	set file( vcf ), file( tbi ) from filtered_snps_first

	output:
	file( "filterd.allBP.vcf.gz" ) into filtered_snps

	script:
	"""
	module load openssl1.0.2

	vcftools \
		--gzvcf ${vcf} \
		--max-missing-count 17 \
		--stdout  \
		--recode | \
		bgzip > filterd.allBP.vcf.gz
	"""
}