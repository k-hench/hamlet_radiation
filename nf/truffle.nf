Channel
	.fromPath("phased_mac2.no_outgroup.vcf.gz")
	.set{ genotypes_ch }

/*
Channel
	.from([[ 50000, 20000, 4 ],
	       [ 25000, 10000, 7 ],
	       [ 15000, 7500, 10 ],
	       [ 10000, 5000, 8 ], 
	       [ 7500, 3000, 9 ],
		   [ 5000, 2000, 3 ]])
	.set{ seq_sizes_ch }
*/

Channel
	.from([[ 25000, 10000, 7 ],
	       [ 15000, 7500, 10 ],
	       [ 10000, 5000, 8 ]])
	.set{ seq_sizes_ch }

Channel
	.from([["direct", ""],
	       ["filter", "LG08"]])
	.set{ filtermode_ch }

process run_truffle {
	publishDir "output/", mode: 'copy'

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
}
