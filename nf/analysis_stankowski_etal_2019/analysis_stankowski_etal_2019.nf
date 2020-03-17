#!/usr/bin/env nextflow
Channel
	.fromFilePairs("../../ressources/other_studies/all_9_taxa_G1_vars.bgziped.vcf.{gz,gz.tbi}")
	.set{ vcf_fl }

Channel
	.fromPath("../../ressources/other_studies/stankowski_etal_2019_spec_short.tsv")
	.set{ pop_fl }

Channel.from( [[1, "ari"], [2, "aur"], [3, "cal"], [4, "cle"], [5, "gra"], [6, "lon"], [7, "par"], [8, "pur"], [9, "puy"]] ).into{ specs_ch1; specs_ch2 }

all_fst_pairs_ch = pop_fl
	.combine(specs_ch1)
	.combine(specs_ch2)
	.filter{ it[1] < it[3] }
	.map{ it[0,2,4]}
	.combine(vcf_fl)

process fst_run {
	label 'L_32g2h_fst_run'
	publishDir "../../ressources/other_studies/stankowski_logs/", mode: 'copy'

	input:
	set file( pop ), val( spec1 ), val( spec2 ), file( vcf_idx ), file( vcf ) from all_fst_pairs_ch

	output:
	file( "${spec1}-${spec2}.log" ) into fst_logs

	script:
	"""
	grep ${spec1} ${pop} > pop1.txt
	grep ${spec2} ${pop} > pop2.txt

	vcftools --gzvcf ${vcf[0]} \
	   --weir-fst-pop pop1.txt \
	   --weir-fst-pop pop2.txt \
	   --out ${spec1}-${spec2} 2> ${spec1}-${spec2}.log
	"""
}

process fst_collect {
	label 'L_loc_fst_collect'
	publishDir "../../ressources/other_studies/", mode: 'copy'

	input:
	set file( log ) from fst_logs.collect()

	output:
	file( "stankowski_etal_2019.tsv" ) into fst_summary

	script:
	"""
	grep "mean Fst estimate" \$BASE_DIR/ressources/other_studies/stankowski_logs/*.log | \
	  sed "s/.log:Weir and Cockerham mean Fst estimate: /\\t/; s=^\$BASE_DIR/ressources/other_studies/stankowski_logs/==" \
		> stankowski_etal_2019.tsv
	"""
}

