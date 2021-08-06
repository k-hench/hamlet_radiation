Channel
	.fromPath("phased_mac2.no_outgroup.vcf.gz")
	.set{ genotypes_ch }

Channel
	.from([[ 50000, 20000, 4 ],
	       [ 25000, 10000, 7 ],
	       [ 15000, 7500, 10 ],
	       [ 10000, 5000, 8 ], 
	       [ 7500, 3000, 9 ],
		   [ 5000, 2000, 3 ]])
	.set{ seq_sizes_ch }


process run_truffle {
	publishDir "output/", mode: 'copy'

	input:
	set file( vcf ), val( sz1 ), val( sz2 ), val( sz3 ) from genotypes_ch.combine( seq_sizes_ch )

	output:
	file( "no_outgr_ld_${sz3}.ibd.tsv" ) into truffle_result

	script:
	"""
	echo "${sz1}_${sz2}"

	truffle \
		--vcf ${vcf} \
		--segments \
		--nofiltering \
		--mindist 5000 \
		--ibs1markers ${sz1} \
		--ibs2markers ${sz2} \
		--out no_outgr_ld_${sz3} \
		--cpu 8
	
	sed 's/^\\s*//g; s/\\s\\+/\\t/g' no_outgr_ld_${sz3}.ibd > no_outgr_ld_${sz3}.ibd.tsv
	"""
}
