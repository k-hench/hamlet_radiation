#!/usr/bin/env nextflow
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

Channel
	.from("windows")
	.set{ window_ch }

Channel
	.from( 50 , 200, 500 )
	.set{ snpnr_ch }

/*mute bed location after test*/
process prep_vcf {
	label "L_20g2h_prep_vcf"
	module "openssl1.0.2"

	input:
	set vcfId, file( vcf ) from vcf_ch

	output:
	set file( "graph.vcf.gz" ), file( "samples_no_out.txt" ) into ( vcf_no_outgr_ch, vcf_no_outgr_slide_ch )

	script:
	"""
	tabix -p vcf -f ${vcf[0]}

	vcfsamplenames ${vcf[0]} | \
		grep -v 'tor\\|tab\\|flo' > samples_no_out.txt

	echo -e "chrom\\tchromStart\\tchromEnd\\nLG01\\t0\\t6472566" > LG01_quarter.bed

	vcftools --gzvcf ${vcf[0]} \
		--bed LG01_quarter.bed \
		--keep samples_no_out.txt \
		--recode \
		--stdout | bgzip >  graph.vcf.gz
	"""
}

snpnr_ch
	.combine( vcf_no_outgr_ch )
	.into{ snpnr_vcf_bg_ch; snpnr_vcf_slide_ch }

process prep_background_graphs {
	label "L_20g2h_bg_graph_prep"
	publishDir "../../2_analysis/popgraphs/data/popgr", mode: 'copy' , pattern: "popgr.background*"
	module "R3.5.2"

	input:
	set val( nsnps ), file( vcf ), file( sample_file ) from snpnr_vcf_bg_ch

	output:
	set val( "bg" ), val( nsnps ), file( "popgr.background*" ), file( sample_file ) into ( bg_graphs_ch )

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/network_generate_background_graphs.R ${nsnps} \
		\$BASE_DIR/R/network_functions.R \
		\$BASE_DIR/R/project_config.R \
		/software/KBIN/downsamplevcf.jar \
		${vcf} \
		${sample_file}
	"""
}

bg_graphs_ch
	.groupTuple()
	.map{[it[0], it[1], it[2], it[3][0]]}
	.set{ bg_graphs_organized_ch }

process summarise_background_graphs {
	label "L_20g2h_bg_graph_summary"
	publishDir "../../2_analysis/popgraphs/figures", mode: 'copy' , pattern: "*.png"
	publishDir "../../2_analysis/popgraphs/data", mode: 'copy' , pattern: "*.tsv.gz"
	module "R3.5.2"

	input:
	set val( bg ), val( n_snps ), file( popgr ), file( sample_file ) from bg_graphs_organized_ch

	output:
	set file( "network_background.png" ), file( "network_background_data.tsv.gz" ) into ( bg_done_ch )

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/network_summarise_background_graphs.R \
		\$BASE_DIR/R/project_config.R \
		\$BASE_DIR/R/network_functions.R \
		\$BASE_DIR/2_analysis/popgraphs/data/popgr/

	gzip network_background_data.tsv
	"""
}


/* un-head prep script */
process prepare_sliding_graphs {
	label "L_20g2h_prep_slide_graphs"

	input:
	set val( n_snps ), file( vcf ), file( sample_file ) from snpnr_vcf_slide_ch

	output:
	file( "*.snp_windows.${n_snps}.tsv" ) into ( slide_windows_ch )
	set val( n_snps ), file( "*.snp_windows.${n_snps}.tsv" ) into ( slide_windows_file )

	script:
	"""
	BASE_NAME=\$(echo ${vcf} | sed 's/.vcf.gz//g')
	\$BASE_DIR/sh/vcf_to_bed_n_snps_windows.sh ${vcf} ${n_snps}
	"""
}

slide_windows_ch
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ lg:row.lg, start:row.start, end:row.end, n_snps:row.n_snps, win_id:row.win_id, gwin_id:row.gwin_id ] }
	.combine( vcf_no_outgr_slide_ch )
	.set { slide_windows_split_ch }

process sliding_graphs {
	label "L_20g15m_slide_graphs"
	publishDir "../../2_analysis/popgraphs/data/popgr", mode: 'copy' , pattern: "popgr.LG*"
	module "R3.5.2"

	input:
	set  val( tab_input ), file( vcf ), file( sample_file ) from slide_windows_split_ch

	output:
	set val( tab_input.n_snps ), file( "popgr.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.rda" ) into ( slide_popgraphs_ch )

	script:
	"""
	echo -e "chrom\\tchromStart\\tchromEnd\\n${tab_input.lg}\\t${tab_input.start}\\t${tab_input.end}" > slide.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.bed

	vcftools --gzvcf ${vcf} \
		--bed slide.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.bed \
		--recode \
		--stdout | gzip > slide.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/network_generate_slide_graphs.R ${tab_input.n_snps} \
		\$BASE_DIR/R/network_functions.R \
		\$BASE_DIR/R/project_config.R \
		${tab_input.lg} \
		${tab_input.start} \
		${tab_input.end} \
		${tab_input.win_id} \
		${tab_input.gwin_id} \
		${sample_file}

	rm slide.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.vcf.gz slide.${tab_input.lg}.${tab_input.n_snps}.${tab_input.win_id}.${tab_input.gwin_id}.bed
	"""
}

slide_popgraphs_ch
	.groupTuple()
	.map{[it[0].toString(), it[1]]}
	.join( slide_windows_file.map{[it[0].toString(), it[1]]})
	.set{ slide_windows_collect }

process sliding_summmary {
	label "L_20g2h_slide_summary"
	publishDir "../../2_analysis/popgraphs/figures", mode: 'copy' , pattern: "*.png"
	publishDir "../../2_analysis/popgraphs/data", mode: 'copy' , pattern: "*.tsv.gz"
	module "R3.5.2"

	input:
	set  val( n_snps ), file( popgr ), file( windows ) from slide_windows_collect

	output:
	set file( "network_slide.${n_snps}.png" ), file( "network_slide_data.${n_snps}.tsv.gz" ) into ( slide_done_ch )

	script:
	"""
	 Rscript --vanilla \$BASE_DIR/R/network_summarise_slide_graphs.R \
		 ${n_snps} \
		 \$BASE_DIR/R/network_functions.R \
		 \$BASE_DIR/R/project_config.R \
		 ${windows}

	gzip network_slide_data.${n_snps}.tsv
	"""
}
