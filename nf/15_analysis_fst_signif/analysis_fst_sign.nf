#!/usr/bin/env nextflow
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations; vcf_genepop_SNP }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

Channel
	.from( "whg", "subset_non_diverged")
	.set{ subset_type_ch }

Channel
	.fromPath( "../../2_analysis/summaries/fst_outliers_998.tsv" )
	.set{ outlier_tab }

// git 3.3
// attach genotypes to location
locations_ch
	.combine( vcf_locations )
	.combine( outlier_tab )
	.combine( subset_type_ch )
	.set{ vcf_location_combo }

process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ), file( outlier_tab ), val( subset_type ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.${subset_type}.vcf.gz" ), file( "${loc}.${subset_type}.vcf.gz.tbi" ), file( "${loc}.${subset_type}.pop" ), val( subset_type ) into ( vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3 )

	script:
	"""
	if [ "${subset_type}" == "subset_non_diverged" ];then
		awk -v OFS="\\t" '{print \$2,\$3,\$4}' ${outlier_tab} > diverged_regions.bed 
		SUBSET="--exclude-bed diverged_regions.bed"
	else
		SUBSET=""
	fi

	vcfsamplenames ${vcf[0]} | \
		grep ${loc} | \
		grep -v tor | \
		grep -v tab > ${loc}.${subset_type}.pop

	vcftools --gzvcf ${vcf[0]} \
		\$SUBSET \
		--keep ${loc}.${subset_type}.pop \
		--mac 3 \
		--recode \
		--stdout | bgzip > ${loc}.${subset_type}.vcf.gz
	
	tabix ${loc}.${subset_type}.vcf.gz
	"""
}

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

process fst_run {
	label 'L_32g1h_fst_run'

	input:
	set val( loc ), file( vcf ), file( vcfidx ), file( pop ), val( subset_type ), val( spec1 ), val( spec2 ) from all_fst_pairs_ch

	output:
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), file( "*_random_fst_a00.tsv" ) into rand_header_ch
	set val( "${spec1}${loc}-${spec2}${loc}_${subset_type}" ), val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( "col1.pop" ), file( "prep.pop" ) into rand_bocy_ch

	script:
	"""
	echo -e "0000\treal_pop" > idx.txt

	vcfsamplenames ${vcf} | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}'  > prep.pop
	grep ${spec1} ${pop} > pop1.txt
	grep ${spec2} ${pop} > pop2.txt
	
	vcftools --gzvcf ${vcf} \
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
	.into{ singles_ch; tens_ch }

singles_ch
	.combine(tens_ch)
	.map{ it[0]+it[1] }
	.toSortedList()
	.flatten()
	.set{ sub_pre_ch }

process random_bodies {
	label 'L_32g6h_fst_run'

	input:
	set val( run ), val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( col1 ), file( prepop ), val( pre ) from rand_bocy_ch.combine(sub_pre_ch)

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
// =======================
// Genepop section
Channel
	.from(["all", "all"], ["hamlets", "all"], ["bel", "loc"], ["hon", "loc"], ["pan", "loc"])
	.set{ genepop_config_ch }

process thin_vcf_genepop {
	label "L_20g2h_subset_vcf"

	input:
	set vcfId, file( vcf ), val( pop ), val( type ) from vcf_genepop_SNP.map{ [it[0].minus(".vcf"), it[1]]}.combine(genepop_config_ch)

	output:
	set val( "${pop}" ) , file( "*_genepop_pops.txt" ) into genepop_prep_ch

	script:
	"""
	module load Java/8.112

	vcfsamplenames ${vcf[0]} | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}' > prep.pop

	if [ "${type}" == "loc" ];then
		grep ${pop} prep.pop | \
		grep -v tor | \
		grep -v tab > pop.txt

		cut -f 1 pop.txt > keepers.txt
		SUBSET="--keep keepers.txt"
	elif [ "${pop}" == "hamlets" ]; then
		grep -v flo prep.pop | \
		grep -v tor | \
		grep -v tab > pop.txt

		cut -f 1 pop.txt > keepers.txt
		SUBSET="--keep keepers.txt"
	else
		mv prep.pop pop.txt
		SUBSET=""
	fi

	# trim linkage last (seprartely, unsure about order within vcftools)
	vcftools \
		--gzvcf ${vcf[0]} \
		--mac 3 \
		\$SUBSET \
		--recode --stdout | \
		vcftools --vcf - \
			--thin 10000 \
			--recode \
			--stdout > ${pop}_sub.vcf

	java -jar \$SFTWR/PGDSpider/PGDSpider2-cli.jar \
		-inputfile ${pop}_sub.vcf \
		-outputfile ${pop}_genepop.txt \
		-spid \$BASE_DIR/ressources/vcf2gp.spid

	sed 's/[A-Za-z0-9_-]*\\([a-z]\\{6\\}\\) ,/\\1 ,/' ${pop}_genepop.txt > ${pop}_genepop_pops.txt
	rm ${pop}_genepop.txt ${pop}_sub.vcf
	"""
}

process run_genepop {
	label "L_32g48h_run_genepop"
	publishDir "../../2_analysis/fst_signif", mode: 'copy' 

	input:
	set val( pop ), file( gp_in ) from genepop_prep_ch

	output:
	set val( pop ), file( "*GE" ), file( "*GE2" ) into genepop_output_ch

	script:
	"""
	Genepop BatchNumber=20 GenepopInputFile=${gp_in} MenuOptions=3.1,3.2 Mode=Batch
	"""
}

process catch_genepop {
	label "L_32g48h_run_genepop"
	publishDir "../../2_analysis/fst_signif", mode: 'copy' 

	input:
	set val( pop ), file( GE ), file( GE2 ) from genepop_output_ch

	output:
	file( "genepop_summary_${pop}.tsv.gz" ) into genepop_catch_ch

	script:
	"""
	tail -n 140 ${GE2} | \
	  grep -v "^--\\|^Normal" | \
	  grep "^[A-Za-z]" | \
	  sed 's/        & /-/; s/y s/y_s/' | \
	  awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' | \
	  gzip > genepop_summary_${pop}.tsv.gz
	"""
}