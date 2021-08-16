#!/usr/bin/env nextflow
// This pipeline includes the analysis calculating the pair-wise IBD

// git 16.1
// load genotypes
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ genotypes_raw_ch }

// git 16.2
// drop outgroups
process drop_outgroups {
	
	input:
	set vcfId, file( vcf ) from genotypes_raw_ch

	output:
	file( "phased_mac2.no_outgroup.vcf.gz" ) into genotypes_ch

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep "tor\\|tab\\|flo" > outrgr.pop

	vcftools --gzvcf ${vcf[0]} \
		--remove outrgr.pop \
		--recode \
		--stdout | bgzip > phased_mac2.no_outgroup.vcf.gz
	"""
}

// git 16.3
// Set IBD fragment sizes
Channel
	.from([[ 25000, 10000, 7 ],
	       [ 15000, 7500, 10 ],
	       [ 10000, 5000, 8 ]])
	.set{ seq_sizes_ch }

// git 16.4
// Set filter mode
Channel
	.from([["direct", ""],
	       ["filter", "LG08"],
	       ["bed", "95"]])
	.set{ filtermode_ch }

// git 16.5
// run truffle
process run_truffle {
	publishDir "../../2_analysis/ibd/", mode: 'copy'

	input:
	set file( vcf ), val( sz1 ), val( sz2 ), val( sz3 ), val( mode ), val( excluding ) from genotypes_ch.combine( seq_sizes_ch ).combine( filtermode_ch )

	output:
	set file( "no_outgr_${mode}${excluding}_${sz3}.ibd.tsv" ), file( "no_outgr_${mode}${excluding}_${sz3}.segments.tsv" ) into truffle_result

	script:
	if( mode == 'direct' )
		"""
		truffle \
			--vcf ${vcf} \
			--segments \
			--nofiltering \
			--ibs1markers ${sz1} \
			--ibs2markers ${sz2} \
			--out no_outgr_${mode}${excluding}_${sz3} \
			--cpu 8
		
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.ibd > no_outgr_${mode}${excluding}_${sz3}.ibd.tsv
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.segments > no_outgr_${mode}${excluding}_${sz3}.segments.tsv
		"""
	else if( mode == 'filter' )
		"""
		vcftools --gzvcf ${vcf} \
			--not-chr ${excluding} \
			--recode \
			--stdout | bgzip > tmp.vcf.gz 
		
		truffle \
			--vcf tmp.vcf.gz  \
			--segments \
			--nofiltering \
			--ibs1markers ${sz1} \
			--ibs2markers ${sz2} \
			--out no_outgr_${mode}${excluding}_${sz3} \
			--cpu 8
		
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.ibd > no_outgr_${mode}${excluding}_${sz3}.ibd.tsv
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.segments > no_outgr_${mode}${excluding}_${sz3}.segments.tsv
		rm tmp.vcf.gz 
		"""
	else if( mode == 'bed' )
		"""
		vcftools --gzvcf ${vcf} \
			--exclude-bed ../../../../../ressources/plugin/idb_above_${excluding}.bed \
			--recode \
			--stdout | bgzip > tmp.vcf.gz 
		
		truffle \
			--vcf tmp.vcf.gz  \
			--segments \
			--nofiltering \
			--ibs1markers ${sz1} \
			--ibs2markers ${sz2} \
			--out no_outgr_${mode}${excluding}_${sz3} \
			--cpu 8
		
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.ibd > no_outgr_${mode}${excluding}_${sz3}.ibd.tsv
		sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_${mode}${excluding}_${sz3}.segments > no_outgr_${mode}${excluding}_${sz3}.segments.tsv
		rm tmp.vcf.gz 
		"""
}

// Duration    : 45m 3s
// CPU hours   : 2.2
// Succeeded   : 3
