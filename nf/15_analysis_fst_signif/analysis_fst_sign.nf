#!/usr/bin/env nextflow
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations; vcf_genepop_SNP }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

// git 3.3
// attach genotypes to location
locations_ch
	.combine( vcf_locations )
	.set{ vcf_location_combo }

process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.vcf.gz.tbi" ), file( "${loc}.pop" ) into ( vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3 )

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

Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

bel_pairs_ch = Channel.from( "bel" )
	.join( vcf_loc_pair1 )
	.combine(bel_spec1_ch)
	.combine(bel_spec2_ch)
	.filter{ it[4] < it[6] }
	.map{ it[0,1,2,3,5,7]}
hon_pairs_ch = Channel.from( "hon" )
	.join( vcf_loc_pair2 )
	.combine(hon_spec1_ch)
	.combine(hon_spec2_ch)
	.filter{ it[4] < it[6] }
	.map{ it[0,1,2,3,5,7]}
pan_pairs_ch = Channel.from( "pan" )
	.join( vcf_loc_pair3 )
	.combine(pan_spec1_ch)
	.combine(pan_spec2_ch)
	.filter{ it[4] < it[6] }
	.map{ it[0,1,2,3,5,7]}
bel_pairs_ch.concat( hon_pairs_ch, pan_pairs_ch  ).set { all_fst_pairs_ch }

process fst_run {
	label 'L_32g1h_fst_run'

	input:
	set val( loc ), file( vcf ), file( vcfidx ), file( pop ), val( spec1 ), val( spec2 ) from all_fst_pairs_ch

	output:
	set val( "${spec1}${loc}-${spec2}${loc}" ), val( loc ), file( "*_random_fst.tsv" ) into rand_header_ch
	set val( "${spec1}${loc}-${spec2}${loc}" ), val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( "col1.pop" ), file( "prep.pop" ) into rand_bocy_ch

	script:
	"""
	vcfsamplenames ${vcf} | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}'  > prep.pop
	grep ${spec1} ${pop} > pop1.txt
	grep ${spec2} ${pop} > pop2.txt
	
	vcftools --gzvcf ${vcf} \
		--weir-fst-pop pop1.txt \
		--weir-fst-pop pop2.txt \
		--stdout 2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	paste idx.txt fst.tsv >> ${spec1}${loc}_${spec2}${loc}_random_fst_00.tsv

	rm fst.tsv fst.log pop1.txt pop2.txt tmp.txt idx.txt

	awk '{print \$1}' prep.pop > col1.pop
	"""
}

Channel
	.from( ('0'..'9'))
	.set{ sub_pre_ch }

rand_bocy_ch.combine(sub_pre_ch).println()
/*	.into{ singles_ch; tens_ch }

singles_ch
	.combine(tens_ch)
	.map{ it[0]+it[1] }
	.toSortedList()
	.flatten()
	.set{ sub_pre_ch }

process split_allBP {
	label 'L_32g30m_fst_run'

	input:
	set val( run ), val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( col1 ), file( prepop ), val( pre ) rand_bocy_ch.combine(sub_pre_ch)

	output:
	file("*txt") into out_ch

	script:
	"""
	for k in {01..99}; do
	echo "Iteration_"\$k
	echo -e "\$k\trandom" > idx.txt

	awk '{print \$2}' ${preppop} | shuf > col2.pop # premutation happens here
	paste ${col1} col2.pop > rand.pop

	grep "${spec1}\$" rand.pop > r_pop1.pop
	grep "${spec2}\$" rand.pop > r_pop2.pop

	vcftools --gzvcf test.vcf.gz \
		--weir-fst-pop r_pop1.pop \
		--weir-fst-pop r_pop2.pop \
		--stdout  2> fst.log 1> tmp.txt

	grep "^Weir" fst.log | sed 's/.* //' | paste - - > fst.tsv
	paste idx.txt fst.tsv >> ${spec1}${loc}_${spec2}${loc}_random_fst_${pre}.tsv

	rm fst.tsv fst.log rand.pop col2.pop r_pop1.pop r_pop2.pop tmp.txt 
	done
	"""
}
*/

// =======================
// Genepop section
process thin_vcf_genepop {
	label "L_20g2h_subset_vcf"

	input:
	set vcfId, file( vcf ) from vcf_genepop_SNP.map{ [it[0].minus(".vcf"), it[1]]}

	output:
	set vcfId, file( "${vcfId}_genepop_pops.txt" ) into genepop_prep_ch

	script:
	"""
	module load Java/8.112

	vcfsamplenames ${vcf[0]} | \
		awk '{print \$1"\\t"substr(\$1, length(\$1)-5, length(\$1))}' > pop.txt

	vcftools \
		--gzvcf ${vcf[0]} \
		--thin 5000 \
		--recode --stdout | \
		vcfrandomsample -r 0.0043 -p 42 > ${vcfId}_sub.vcf # 0.43% ~ 100k SNPs

	java -jar \$SFTWR/PGDSpider/PGDSpider2-cli.jar \
		-inputfile ${vcfId}_sub.vcf \
		-outputfile ${vcfId}_genepop.txt \
		-spid \$BASE_DIR/ressources/vcf2gp.spid

	sed 's/[A-Za-z0-9_-]*\\([a-z]\\{6\\}\\) ,/\\1 ,/' ${vcfId}_genepop.txt > ${vcfId}_genepop_pops.txt
	"""
}
/*
process run_genepop {
	label "L_20g2h_run_genepop"

	input:
	set vcfId, file( gp_in ) into genepop_prep_ch

	output:
	set file( "*GE" ) into genepop_prep_ch

	script:
	"""
	Genepop BatchNumber=20 GenepopInputFile=${gp_in} MenuOptions=3.1,3.2 Mode=Batch
	"""
}
*/
/*

# module load Java/8.112
java -jar /software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile test.vcf -outputfile test_genepop_geno.txt -spid vcf2gp_no_pops.spid

sed 's/[A-Za-z0-9_-]*\([a-z]\{6\}\) ,/\1 ,/' test_genepop_geno.txt > test_genepop_geno_pops.txt

Genepop BatchNumber=20 GenepopInputFile=test_genepop_geno_pops.txt MenuOptions=3.1,3.2 Mode=Batch

*/