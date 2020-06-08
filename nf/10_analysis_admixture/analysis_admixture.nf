#!/usr/bin/env nextflow
// git 10.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

// git 10.2
// Set different k values for the admixture analysis
Channel
	.from( 2..15 )
	.set{ k_ch }

// git 10.3
// load Fst outlier regions
Channel
	.fromPath("../../ressources/plugin/poptrees/outlier.bed")
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ chrom:row.chrom, start:row.start, end:row.end, gid:row.gid ] }
	.combine( vcf_ch )
	.set{ vcf_admx }

// git 10.4
// subset genotypes to the outlier region and reformat
process plink12 {
	label 'L_20g2h_plink12'
	tag "${grouping.gid}"

	input:
	set val( grouping ),  val( vcfidx ), file( vcf ) from vcf_admx

	output:
	set val( grouping ), file( "hapmap.*.ped" ), file( "hapmap.*.map" ), file( "hapmap.*.nosex" ), file( "pop.txt" ) into admx_plink

	script:
	"""
	echo -e "CHROM\\tSTART\\tEND" > outl.bed
	echo -e "${grouping.chrom}\\t${grouping.start}\\t${grouping.end}" >> outl.bed

	vcfsamplenames ${vcf[0]} | \
		grep -v "tor\\|tab\\|flo" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > pop.txt

	vcftools \
		--gzvcf ${vcf[0]} \
		--keep pop.txt \
		--bed outl.bed \
		--plink \
		--out admx_plink

	plink \
		--file admx_plink \
		--recode12 \
		--out hapmap.${grouping.gid}
	"""
}

// git 10.5
// combine genoutype subsets with k values
admx_prep  = k_ch.combine( admx_plink )

// git 10.6
// run admixture
process admixture_all {
	label 'L_20g4h_admixture_all'
	publishDir "../../2_analysis/admixture/", mode: 'copy'
	tag "${grouping.gid}.${k}"

	input:
	set  val( k ), val( grouping ), file( ped ), file( map ), file( nosex ), file( pop ) from admx_prep

	output:
	set val( "dummy" ), file( "*.out" ), file( "*.Q" ), file( "*.txt" ) into admx_log

	script:
	"""
	mv ${pop} pop.${grouping.gid}.${k}.txt
	admixture --cv ${ped} ${k} | tee log.${grouping.gid}.${k}.out
	"""
}
