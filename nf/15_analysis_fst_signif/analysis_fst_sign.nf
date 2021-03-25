#!/usr/bin/env nextflow
// git 15.1
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations; vcf_adapt }

// git 15.2
Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

// git 15.3
Channel
	.from( "whg", "subset_non_diverged")
	.into{ subset_type_ch; subset_type_ch2 }

// git 15.4
Channel
	.fromPath( "../../2_analysis/summaries/fst_outliers_998.tsv" )
	.into{ outlier_tab; outlier_tab2 }

// git 15.5
// attach genotypes to location
locations_ch
	.combine( vcf_locations )
	.combine( outlier_tab )
	.set{ vcf_location_combo }

// git 15.6
process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ), file( outlier_tab ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.vcf.gz.tbi" ), file( "${loc}.pop" ), file( outlier_tab ) into ( vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3 )

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
		--stdout | bgzip > ${loc}.vcf.gz
	
	tabix ${loc}.vcf.gz
	"""
}

// git 15.7
Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

bel_pairs_ch = Channel.from( "bel" )
	.join( vcf_loc_pair1 )
	.combine(bel_spec1_ch)
	.combine(bel_spec2_ch)
	.filter{ it[5] < it[7] }
	.map{ it[0,1,2,3,4,6,8]}
hon_pairs_ch = Channel.from( "hon" )
	.join( vcf_loc_pair2 )
	.combine(hon_spec1_ch)
	.combine(hon_spec2_ch)
	.filter{ it[5] < it[7] }
	.map{ it[0,1,2,3,4,6,8]}
pan_pairs_ch = Channel.from( "pan" )
	.join( vcf_loc_pair3 )
	.combine(pan_spec1_ch)
	.combine(pan_spec2_ch)
	.filter{ it[5] < it[7] }
	.map{ it[0,1,2,3,4,6,8]}
bel_pairs_ch.concat( hon_pairs_ch, pan_pairs_ch  ).set { all_fst_pairs_ch }

// git 15.8
process fst_run {
	label 'L_32g1h_fst_run'

	input:
	set val( loc ), file( vcf ), file( vcfidx ), file( pop ), file( outlier_tab ), val( spec1 ), val( spec2 ), val( subset_type ) from all_fst_pairs_ch.combine( subset_type_ch )

	output:
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), file( "*_random_fst_a00.tsv" ) into rand_header_ch
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), val( loc ), val( spec1 ), val( spec2 ), file( "${loc}.${subset_type}.vcf.gz" ), file( "col1.pop" ), file( "prep.pop" ) into rand_body_ch

	script:
	"""
	if [ "${subset_type}" == "subset_non_diverged" ];then
		awk -v OFS="\\t" '{print \$2,\$3,\$4}' ${outlier_tab} > diverged_regions.bed 
		SUBSET="--exclude-bed diverged_regions.bed"
	else
		SUBSET=""
	fi

	vcftools --gzvcf ${vcf[0]} \
		\$SUBSET \
		--recode \
		--stdout | bgzip > ${loc}.${subset_type}.vcf.gz

	tabix ${loc}.${subset_type}.vcf.gz

	echo -e "0000\treal_pop" > idx.txt

	vcfsamplenames ${loc}.${subset_type}.vcf.gz | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}'  > prep.pop
	grep ${spec1} ${pop} > pop1.txt
	grep ${spec2} ${pop} > pop2.txt
	
	vcftools --gzvcf ${loc}.${subset_type}.vcf.gz \
		--weir-fst-pop pop1.txt \
		--weir-fst-pop pop2.txt \
		--stdout 2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	echo -e "idx\\ttype\\tmean_fst\\tweighted_fst" > ${spec1}${loc}-${spec2}${loc}_${subset_type}_random_fst_a00.tsv
	paste idx.txt fst.tsv >> ${spec1}${loc}-${spec2}${loc}_${subset_type}_random_fst_a00.tsv

	rm fst.tsv fst.log pop1.txt pop2.txt tmp.txt idx.txt

	awk '{print \$1}' prep.pop > col1.pop
	"""
}


Channel
	.from( ('0'..'9'))
	.map{ "0" + it }.into{ sub_pre_ch; sub_pre_ch2 }
/*
	.into{ singles_ch; tens_ch }

// git 15.9
singles_ch
	.combine(tens_ch)
	.map{ it[0]+it[1] }
	.toSortedList()
	.flatten()
	.into{ sub_pre_ch; sub_pre_ch2 }*/

// git 15.10
process random_bodies {
	label 'L_32g6h_fst_run'

	input:
	set val( run ), val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( col1 ), file( prepop ), val( pre ) from rand_body_ch.combine(sub_pre_ch)

	output:
	set val( run ), file("*_random_fst_b${pre}.tsv") into rand_body_out_ch

	script:
	"""
	for k in {00..99}; do
	echo "Iteration_"\$k
	echo -e "${prepop}\$k\trandom" > idx.txt

	awk '{print \$2}' ${prepop} | shuf > col2.pop # premutation happens here
	paste ${col1} col2.pop > rand.pop

	grep "${spec1}${loc}\$" rand.pop > r_pop1.pop
	grep "${spec2}${loc}\$" rand.pop > r_pop2.pop

	vcftools --gzvcf ${vcf} \
		--weir-fst-pop r_pop1.pop \
		--weir-fst-pop r_pop2.pop \
		--stdout  2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	paste idx.txt fst.tsv >> ${run}_random_fst_b${pre}.tsv

	rm fst.tsv fst.log rand.pop col2.pop r_pop1.pop r_pop2.pop tmp.txt 
	done
	"""
}

// git 15.11
process compile_random_results {
	label 'L_20g2h_compile_rand'
	publishDir "../../2_analysis/fst_signif/random", mode: 'copy' 

	input:
	set val( run ), file( body ), file( head ) from rand_body_out_ch.groupTuple().join(rand_header_ch, remainder: true)

	output:
	file("${run}_random_fst.tsv.gz") into random_lists_result

	script:
	"""
	cat ${head} > ${run}_random_fst.tsv
	cat ${body} >> ${run}_random_fst.tsv
	gzip ${run}_random_fst.tsv
	"""
}

// same for adaptation -----------
// git 15.12
Channel
	.from( "nig", "pue", "uni")
	.set{ species_ch }

// git 15.13
// define location set
Channel.from( [[1, "bel"], [2, "hon"], [3, "pan"]]).into{ locations_ch_1;locations_ch_2 }

// git 15.14
// create location pairs
locations_ch_1
	.combine(locations_ch_2)
	.filter{ it[0] < it[2] }
	.map{ it[1,3]}
	.combine( species_ch )
	.combine( vcf_adapt )
	.combine( outlier_tab2 )
	.combine( subset_type_ch2 )
	.set{ vcf_location_combo_adapt }

// git 15.15
process fst_run_adapt {
	label 'L_32g1h_fst_run'

	input:
	set val( loc1 ), val( loc2 ), val( spec ), val( vcf_indx) , file( vcf ), file( outlier_tab ), val( subset_type ) from vcf_location_combo_adapt

	output:
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), file( "*_random_fst_a00.tsv" ) into rand_header_adapt_ch
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), val( spec ), val( loc1 ), val( loc2 ), file( "${loc}.${subset_type}.vcf.gz" ), file( "col1.pop" ), file( "prep.pop" ) into rand_body_adapt_ch

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep ${spec} > ${spec}.pop

	if [ "${subset_type}" == "subset_non_diverged" ];then
		awk -v OFS="\\t" '{print \$2,\$3,\$4}' ${outlier_tab} > diverged_regions.bed 
		SUBSET="--exclude-bed diverged_regions.bed"
	else
		SUBSET=""
	fi

	vcftools --gzvcf ${vcf[0]} \
		\$SUBSET \
		--keep ${spec}.pop \
		--mac 3 \
		--recode \
		--stdout | bgzip > ${spec}.${subset_type}.vcf.gz

	tabix ${spec}.${subset_type}.vcf.gz

	echo -e "0000\treal_pop" > idx.txt

	vcfsamplenames ${spec}.${subset_type}.vcf.gz | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}'  > prep.pop
	grep ${loc1} ${spec}.pop > pop1.txt
	grep ${loc2} ${spec}.pop > pop2.txt
	
	vcftools --gzvcf ${spec}.${subset_type}.vcf.gz \
		--weir-fst-pop pop1.txt \
		--weir-fst-pop pop2.txt \
		--stdout 2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	echo -e "idx\\ttype\\tmean_fst\\tweighted_fst" > ${spec}${loc1}-${spec}${loc2}_${subset_type}_random_fst_a00.tsv
	paste idx.txt fst.tsv >> ${spec}${loc1}-${spec}${loc2}_${subset_type}_random_fst_a00.tsv

	rm fst.tsv fst.log pop1.txt pop2.txt tmp.txt idx.txt

	awk '{print \$1}' prep.pop > col1.pop
	"""
}

// git 15.16
process random_bodies_adapt {
	label 'L_32g6h_fst_run'

	input:
	set val( run ), val( spec ), val( loc1 ), val( loc2 ), file( vcf ), file( col1 ), file( prepop ), val( pre ) from rand_body_adapt_ch.combine(sub_pre_ch2)

	output:
	set val( run ), file("*_random_fst_b${pre}.tsv") into rand_body_out_adapt_ch

	script:
	"""
	for k in {00..99}; do
	echo "Iteration_"\$k
	echo -e "${prepop}\$k\trandom" > idx.txt

	awk '{print \$2}' ${prepop} | shuf > col2.pop # premutation happens here
	paste ${col1} col2.pop > rand.pop

	grep "${spec}${loc1}\$" rand.pop > r_pop1.pop
	grep "${spec}${loc2}\$" rand.pop > r_pop2.pop

	vcftools --gzvcf ${vcf} \
		--weir-fst-pop r_pop1.pop \
		--weir-fst-pop r_pop2.pop \
		--stdout  2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	paste idx.txt fst.tsv >> ${run}_random_fst_b${pre}.tsv

	rm fst.tsv fst.log rand.pop col2.pop r_pop1.pop r_pop2.pop tmp.txt 
	done
	"""
}

// git 15.17
process compile_random_results_adapt {
	label 'L_20g2h_compile_rand'
	publishDir "../../2_analysis/fst_signif/random/adapt", mode: 'copy' 

	input:
	set val( run ), file( body ), file( head ) from rand_body_out_adapt_ch.groupTuple().join(rand_header_adapt_ch, remainder: true)

	output:
	file("${run}_random_fst.tsv.gz") into random_lists_adapt_result

	script:
	"""
	cat ${head} > ${run}_random_fst.tsv
	cat ${body} >> ${run}_random_fst.tsv
	gzip ${run}_random_fst.tsv
	"""
}