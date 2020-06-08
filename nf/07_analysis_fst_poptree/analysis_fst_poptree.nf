#!/usr/bin/env nextflow
// git 7.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_fst }

// git 7.2
// load all possible population pairs
Channel
	.fromPath("../../ressources/plugin/poptrees/all_crosses.tsv")
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ pop1:row.pop1, pop2:row.pop2 ] }
	.set{ crosses_ch }

// git 7.3
// open the focal Fst outlier regions
Channel
	.fromPath("../../ressources/plugin/poptrees/outlier.bed")
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ chrom:row.chrom, start:row.start, end:row.end, gid:row.gid ] }
	.combine( vcf_fst )
	.combine( crosses_ch )
	.set{ crosses_vcf }

// git 7.4
// compute the average Fst for all possible pair within the outlier region
process outlier_fst {
		label "L_loc_collect_fst"
		publishDir "../../2_analysis/fst/poptree/single", mode: 'copy'

		input:
		set val( grouping ),  val( vcfidx ), file( vcf ), val( cross_pop ) from crosses_vcf

		output:
		set val( grouping.gid ), val( cross_pop ), file( "*.fst.tsv" ) into outlier_fst_gid_ch

		script:
		"""
		echo -e "CHROM\\tSTART\\tEND" > outl.bed
		echo -e "${grouping.chrom}\\t${grouping.start}\\t${grouping.end}" >> outl.bed

		vcfsamplenames ${vcf[0]} | \
			grep "${cross_pop.pop1}" > pop1.pop

			vcfsamplenames ${vcf[0]} | \
				grep "${cross_pop.pop2}" > pop2.pop

		vcftools --gzvcf ${vcf[0]} \
			--bed outl.bed \
			--keep pop1.pop \
			--keep pop2.pop \
			--weir-fst-pop pop1.pop \
			--weir-fst-pop pop2.pop \
			--stdout 2> ${cross_pop.pop1}-${cross_pop.pop2}.50k.log | \
			gzip > ${cross_pop.pop1}-${cross_pop.pop2}.fst.tsv.gz

		mFST=\$(grep "Weir and Cockerham mean Fst estimate:" ${cross_pop.pop1}-${cross_pop.pop2}.50k.log | sed 's/Weir and Cockerham mean Fst estimate: //')
		wFST=\$(grep "Weir and Cockerham weighted Fst estimate:" ${cross_pop.pop1}-${cross_pop.pop2}.50k.log | sed 's/Weir and Cockerham weighted Fst estimate: //')

		echo -e "${cross_pop.pop1}-${cross_pop.pop2}\\t\$mFST\\t\$wFST" > ${cross_pop.pop1}-${cross_pop.pop2}.${grouping.gid}.fst.tsv
		"""
	}

// git 7.5
// collect all population pairs within each region and compile Fst table
process outlier_fst_collect {
		label "L_20g2h_outlier_fst"
		publishDir "../../2_analysis/fst/poptree/summary", mode: 'copy'

		input:
		set val( gid ), val( cross_pop ), file( fst ) from outlier_fst_gid_ch.groupTuple()

		output:
		file( "${gid}.fst.all.tsv" ) into outlier_fst_collect_ch

		script:
		"""
		echo -e "run\\tmean_fst\\tweighted_fst" > ${gid}.fst.all.tsv
		cat *.fst.tsv >> ${gid}.fst.all.tsv
		"""
	}

// git 7.6
// neighbour joining happens within the R script R/fig/plot_F4.R
