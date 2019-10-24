#!/usr/bin/env nextflow
// This pipelie includes the recombination anlysis

// git 6.1
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

// git 6.2
Channel
	.from( 1..24 )
	.map{ it.toString().padLeft(2, "0") }
	.set{ lg_ch }

// git 6.3
// split the genotypes by LG and reformat the genotypes
process split_allBP {
	label 'L_20g2h_split_by_lg'
	tag "LG${lg}"

	input:
	set val( lg ), vcfId, file( vcf ) from lg_ch.combine( vcf_ch )

	output:
	set val( lg ), file( "phased_mac2.LG${lg}.vcf.gz" ) into vcf_by_lg_ch

	script:
	"""
	module load openssl1.0.2

	vcftools --gzvcf ${vcf[0]} \
		--chr LG${lg} \
		--recode \
		--stdout | bgzip  > phased_mac2.LG${lg}.vcf.gz
	"""
}

// git 6.4
process fasteprr_s1 {
	label 'L_20g2h_fasteprr_s1'
	tag "LG${lg}"
	module "R3.5.2"

	input:
	set val( lg ), file( vcf ) from vcf_by_lg_ch

	output:
	file( "step1_LG${lg}" ) into step_1_out_ch

	script:
	"""
	mkdir step1_LG${lg}
	Rscript --vanilla \$BASE_DIR/R/fasteprr_step1.R ./${vcf} step1_LG${lg} LG${lg} 50
	"""
}

// git 6.5
process fasteprr_s1_summary {
	label 'L_loc_fasteprr_s1_summmary'

	input:
	file( step1 ) from step_1_out_ch.collect()

	output:
	file( "step1" ) into ( step1_ch1, step1_ch2 )

	script:
	"""
	mkdir step1
	cp step1_LG*/* step1/
	"""
}

// git 6.6
Channel
	.from( 1..250 )
	.map{ it.toString().padLeft(3, "0") }
	.combine( step1_ch1 )
	.set{ step_2_run_ch }

// git 6.7
process fasteprr_s2 {
	label 'L_long_loc_fasteprr_s2'
	tag "run_${idx}"
	module "R3.5.2"

	input:
	set val( idx ), file( step1 ) from step_2_run_ch

	output:
	set val( idx ), file( "step2_run${idx}" ) into step_2_out_ch

	script:
	"""
	mkdir -p step2_run${idx}
	Rscript --vanilla \$BASE_DIR/R/fasteprr_step2.R ${step1} step2_run${idx} ${idx}
	"""
}

step_2_out_ch.into{ step_2_indxs; step_2_files }

step_2_indxs.map{ it[0] }.collect().println()
step_2_files.map{ it[1] }.collect().println()
/*
// git 6.8
process fasteprr_s2_summary {
	label 'L_loc_fasteprr_s2_summmary'

	input:
	set val( idx ), file( step2 ) from step_2_out_ch.collect()

	output:
	file( "step2" ) into ( step2_ch )

	script:
	"""
	mkdir step2

	for k in ${idx};do
		cp -r step2_run\$k/\$k step2/
	done
	"""
}

// git 6.9
process fasteprr_s3 {
	label 'L_32g4h_fasteprr_s3'
	module "R3.5.2"

	input:
	set file( step1 ), file( step2 ) from step1_ch2.combine( step2_ch )

	output:
	file( "step3" ) into step_3_out_ch

	script:
	"""
	mkdir step3
	Rscript --vanilla \$BASE_DIR/R/fasteprr_step3.R ${step1} ${step2} step3
	"""
}

// git 6.10
process fasteprr_s3_summary {
	label 'L_loc_fasteprr_s3_summmary'
	publishDir "../../2_analysis/fasteprr", mode: 'copy'

	input:
	file( step3 ) from step_3_out_ch

	output:
	file( "step4/fasteprr.all.rho.txt.gz" ) into ( step3_ch )

	script:
	"""
	mkdir step4

	# ------ rewriting the fasteprr output into tidy format --------

	for k in  {01..24};do
	    j="LG"\$k;
	    echo \$j;
	    \$BASE_DIR/sh/fasteprr_trans.sh step3/chr_\$j \$j step4/fasteprr.\$j
	done

	# --------- combining all LGs into a single data set -----------
	cd step4

	head -n 1 fasteprr.LG01.rho.txt > fasteprr.all.rho.txt

	for k in {01..24}; do
	echo "LG"\$k
	awk 'NR>1{print \$0}' fasteprr.LG\$k.rho.txt >> fasteprr.all.rho.txt
	done
	gzip fasteprr.all.rho.txt
	cd ..
	"""
}
*/
