process twisst_sample_grouping {
	label "L_loc_msmc_grouping"
	publishDir "../../2_analysis/twisst_crossloc/setup", mode: 'copy'

	output:
	file( "twisst_runs_config.tsv" ) into twisst_table_ch

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/twisst_crossloc_grouping.R
	"""
}

twisst_table_ch
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ runnr:row.runnr, pop1:row.pop1, pop2:row.pop2, pop3:row.pop3, pop4:row.pop4 ] }
	.set { twisst_grouping_ch }

Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

	Channel
		.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
		.map{ "LG" + it }
		.set{ lg_ch }

Channel
	.from(50, 200)
	.combine( twisst_grouping_ch )
	.combine( lg_ch )
	.combine( vcf_ch )
	.set{ twisst_runs_ch }

process subset_vcf_by_location_whg {
		label "L_30g2h4x_subset_vcf_whg"
		publishDir "../../2_analysis/twisst_crossloc/run_${grouping.runnr}/positions/", mode: 'copy', patttern: "*.data.tsv"

		input:
		set val( nSNP ), val( grouping ), val( lg ), val( vcfidx ), file( vcf ) from twisst_runs_ch

		output:
		set val( nSNP ), val( grouping ), val( lg ), file( "*.trees.gz" ), file( "*.data.tsv" ), file("pop_all.*") into twisst_phylos_ch

		script:
		"""
		vcfsamplenames ${vcf[0]} > pop_prep.txt

		grep  ${grouping.pop1} pop_prep.txt > pop1.pop
		grep  ${grouping.pop2} pop_prep.txt > pop2.pop
		grep  ${grouping.pop3} pop_prep.txt > pop3.pop
		grep  ${grouping.pop4} pop_prep.txt > pop4.pop

		cat *.pop | \
			awk -v OTF="\\t" '{print \$1,substr( \$1, length(\$1) - 5, length(\$1) )}' \
			> pop_all.${grouping.runnr}.pop

		vcftools --gzvcf  ${vcf[0]} \
			--keep pop1.pop \
			--keep pop2.pop \
			--keep pop3.pop \
			--keep pop4.pop \
			--chr ${lg} \
			--recode --stdout | \
			gzip > twisst_crossloc.${grouping.runnr}.${lg}.${nSNP}SNPS.vcf.gz

		python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
			-i  twisst_crossloc.${grouping.runnr}.${lg}.${nSNP}SNPS.vcf.gz |
			python \$SFTWR/genomics_general/phylo/phyml_sliding_windows.py \
			--threads 2 \
			--prefix twisst_crossloc.${grouping.runnr}.${lg}.${nSNP}SNPS \
			--windType sites \
			-w ${nSNP} \
			--model GTR \
			--optimise n
		"""
	}


process run_twisst {
	label "L_20g2h_run_twisst"
	publishDir "../../2_analysis/twisst_crossloc/run_${grouping.runnr}/weights/", mode: 'copy'

	input:
	set val( nSNP ), val( grouping ), val( lg ), file( trees ), file( data ), file( groups ) from twisst_phylos_ch

	output:
	file( "*.weights.tsv" ) into final_result

	script:
	"""
	awk -v OFT="\\t" '{print \$1"_A",\$2"\\n"\$1"_B",\$2}' ${groups} \
		> twisst.${grouping.runnr}.pop

	python \$SFTWR/twisst/twisst.py \
		-t ${trees} \
		-w twisst_crossloc.${grouping.runnr}.${lg}.${nSNP}SNPS.weights.tsv \
		-g ${grouping.pop1} \
		-g ${grouping.pop2} \
		-g ${grouping.pop3} \
		-g ${grouping.pop4} \
		--groupsFile ${groups}
	"""
}