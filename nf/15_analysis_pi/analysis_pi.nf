#!/usr/bin/env nextflow
// This pipeline includes the analysis run on the
//   all callable sites data sheet (pi).

// git 15.1
// load genotypes
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.vcf.{gz,gz.tbi}")
	.set{ vcf_pi_ch }

// git 15.2
// initialize LGs
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ lg_pi_ch }

// git 15.3
// set outlier window mode
Channel
	.from( "", "_no_outlier" )
	.set{ outlier_mode_ch }

// git 15.4
// init slining window resolutions
Channel
	.from( 1, 5 )
	.set{ kb_ch }

// git 15.5
// init all sampled populations (for pi)
Channel
	.from('indbel', 'maybel', 'nigbel', 'puebel', 'unibel', 'abehon', 'gumhon', 'nighon', 'puehon', 'ranhon', 'unihon', 'nigpan', 'puepan', 'unipan')
	.combine( vcf_pi_ch )
	.combine( kb_ch )
	.combine( lg_pi_ch )
	.combine( outlier_mode_ch )
	.set{ input_ch }

// git 15.5
// filter genotypes and convert formats
process recode_genotypes {
	label 'L_32g15h_recode'
	tag "${spec}"

	input:
	set val( spec ), vcfId, file( vcf ), val( kb ), val( lg ), val( outlier ) from input_ch

	output:
	set val( spec ), val( kb ), val( lg ), val( outlier ), file( "*.geno.gz" ), file( "pop.txt" ) into geno_ch

	script:
	"""
	if [ "${outlier}" == "_no_outlier" ];then
		tail -n +2 \$BASE_DIR/2_analysis/summaries/fst_outliers_998.tsv | \
			cut -f 2,3,4 > outlier.bed 
		SUBSET="--exclude-bed outlier.bed"
	else
		SUBSET=""
	fi

	vcfsamplenames ${vcf[0]} | \
		grep ${spec} > pop.txt

	vcftools --gzvcf ${vcf[0]} \
		--keep pop.txt \
		--chr LG${lg} \
		\$SUBSET \
		--recode \
		--stdout | bgzip > ${spec}${outlier}.LG${lg}.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i ${spec}${outlier}.LG${lg}.vcf.gz | \
		gzip > ${spec}${outlier}.LG${lg}.geno.gz
	"""
}

// git 15.6
// calculate pi per species
process pi_per_spec {
	label 'L_32g15h_pi'
	tag "${spec}"

	input:
	set val( spec ), val( kb ), val( lg ), val( outlier ), file( geno ), file( pop ) from geno_ch

	output:
	set val( "${spec}${outlier}.${kb}" ), val( kb ), file( "*0kb.csv.gz" ) into pi_lg_ch

	script:
	"""
	awk '{print \$1"\\t"substr(\$1,length(\$1)-5,length(\$1)-1)}' ${pop} > ${spec}.pop

	python \$SFTWR/genomics_general/popgenWindows.py \
		-w ${kb}0000 \
		-s ${kb}000 \
		--popsFile ${spec}.pop \
		-p ${spec} \
		-g ${geno} \
		-o pi.${spec}${outlier}.LG${lg}.${kb}0kb.csv.gz \
		-f phased \
		--writeFailedWindows \
		-T 1
	"""
}

process merge_pi {
	label 'L_32g4h_pi_merge'
	tag "${spec}"
	publishDir "../../2_analysis/pi/${kb[0]}0k", mode: 'copy'

	input:
	set val( spec_outlier_kb ), val( kb ), file( pi ) from pi_lg_ch.groupTuple()

	output:
	file( "pi.${spec_outlier_kb}0kb.tsv.gz" ) into pi_output_ch

	script:
	"""
	echo -e "CHROM\\tBIN_START\\tBIN_END\\tBIN_MID_SITES\\tN_SITES\\tPI" > pi.${spec_outlier_kb}0kb.tsv
	zcat ${pi} | \
		tail -n +2 | \
		grep -v "==>" | \
		grep -v "^\$" | \
		sed -s "s/,/\\t/g" >> pi.${spec_outlier_kb}0kb.tsv
	
	gzip pi.${spec_outlier_kb}0kb.tsv
	"""
}