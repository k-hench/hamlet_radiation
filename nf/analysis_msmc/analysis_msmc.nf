#!/usr/bin/env nextflow
/* create channel of linkage groups */
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.map{ "LG" + it }
	.into{ lg_ch1; lg_ch2; lg_ch3 }

Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_msmc }

/* 1) Msmc section ============== */

Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd_bi-allelic.vcf.{gz,gz.tbi}")
	.set{ vcf_depth }

/* gather depth per individual ----------------------------- */
process gather_depth {
	label 'L_20g2h_split_by_sample'
	publishDir "metadata", mode: 'copy'

	input:
	set vcfID, file( vcf ) from vcf_depth

	output:
	file( "depth_by_sample.txt" ) into depth_ch
	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		 grep -v "tor\\|tab\\|flo" > pop.txt

	vcftools \
		--gzvcf ${vcf[0]} \
		--keep pop.txt \
		--depth \
		--stdout > depth_by_sample.txt
	"""
}

depth_ch
/* create channel out of sequencing depth table */
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ id:row.INDV, sites:row.N_SITES, depth:row.MEAN_DEPTH] }
	.map{ [it.id, it.sites, it.depth] }
	.set { depth_by_sample_ch }

/* create channel from bam files and add sample id */
Channel
	.fromPath( '../../1_genotyping/0_dedup_bams/*.bam' )
	.map{ file ->
				def key = file.name.toString().tokenize('.').get(0)
				return tuple(key, file)}
				.set{ sample_bams }

/* combine sample bams and sequencing depth */
sample_bams
	.join( depth_by_sample_ch )
	.set{ sample_bam_and_depth }

/* multiply the sample channel by the linkage groups */
sample_bam_and_depth
	.combine( vcf_msmc )
	.combine( lg_ch1 )
	.set{ samples_msmc }

/* split vcf by individual ----------------------------- */
process split_vcf_by_individual {
	label 'L_20g15m_split_by_vcf'

	input:
	set val( id ), file( bam ), val( sites ), val( depth ), val( vcf_id ), file( vcf ), val( lg ) from samples_msmc

	output:
	set val( id ), val( lg ), file( bam ), val( depth ), file( "phased_mac2.${id}.${lg}.vcf.gz" ) into ( sample_vcf, sample_vcf2 )

	script:
	"""
	gatk --java-options "-Xmx10G" \
		SelectVariants \
		-R \$REF_GENOME \
		-V ${vcf[0]} \
		-sn ${id} \
		-L ${lg}\
		-O phased_mac2.${id}.${lg}.vcf.gz
	"""
}

process bam_caller {
	label 'L_36g47h_bam_caller'
	publishDir "../../ressources/coverage_masks", mode: 'copy' , pattern: "*.coverage_mask.bed.gz"
	conda "$HOME/miniconda2/envs/py3"

	input:
	set val( id ), val( lg ), file( bam ), val( depth ), file( vcf ) from sample_vcf

	output:
	set val( id ), val( lg ), file( "*.bam_caller.vcf.gz" ), file( "*.coverage_mask.bed.gz" ) into coverage_by_sample_lg

	script:
	"""
	module load openssl1.0.2

	samtools index ${bam}

	samtools mpileup -q 25 -Q 20 -C 50 -u -r ${lg} -f \$REF_GENOME ${bam} | \
		bcftools call -c -V indels | \
		\$BASE_DIR/py/bamHamletCaller.py ${depth} ${id}.${lg}.coverage_mask.bed.gz | \
		gzip -c > ${id}.${lg}.bam_caller.vcf.gz
	"""
}

process generate_segsites {
	label "L_20g15m_msmc_generate_segsites"
	publishDir "../../2_analysis/msmc/segsites", mode: 'copy' , pattern: "*.segsites.vcf.gz"

	input:
	set val( id ), val( lg ), file( bam ), val( depth ), file( vcf ) from sample_vcf2

	output:
	set val( id ), val( lg ), file( "*.segsites.vcf.gz" ), file( "*.covered_sites.bed.txt.gz" ) into segsites_by_sample_lg

	script:
	"""
	zcat ${vcf} | \
		vcfAllSiteParser.py ${id} ${id}.${lg}.covered_sites.bed.txt.gz | \
		gzip -c > ${id}.${lg}.segsites.vcf.gz
	"""
}

process msmc_sample_grouping {
	label "L_loc_msmc_grouping"
	publishDir "../../2_analysis/msmc/setup", mode: 'copy'
	module "R3.5.2"

	output:
	file( "msmc_grouping.txt" ) into msmc_grouping
	file( "msmc_cc_grouping.txt" ) into  cc_grouping

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/sample_assignment_msmc.R \
			\$BASE_DIR/R/distribute_samples_msmc_and_cc.R \
			\$BASE_DIR/R/cross_cc.R \
			\$BASE_DIR/metadata/sample_info.txt \
			msmc
	"""
}

msmc_grouping
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ run:row.msmc_run, spec:row.spec, geo:row.geo, group_nr:row.group_nr, group_size:row.group_size, samples:row.samples ] }
	.set { msmc_runs }
/* wait for bam_caller and generate_segsites to finish: */
/*this '.collect' is only meant to wait until the channel is done,
  files are being redirected via publishDir*/
coverage_by_sample_lg.collect().map{ "coverage done!" }.into{ coverage_done; coverage_cc }
segsites_by_sample_lg.collect().map{ "segsites done!" }.into{ segsites_done; segsites_cc }

lg_ch2
	.combine( msmc_runs )
	.combine( coverage_done )
	.combine( segsites_done )
	.map{[it[0], it[1].run, it[1]]}
	.set{ msmc_grouping_after_segsites }

/* generating MSMC input files (4 or 3 inds per species) ----------- */

process generate_multihetsep {
	label "L_120g40h_msmc_generate_multihetsep"
	publishDir "../../2_analysis/msmc/input/run_${run}", mode: 'copy' , pattern: "*.multihetsep.txt"
	conda "$HOME/miniconda2/envs/py3"

	input:
	/* content msmc_gr: val( msmc_run ), val( spec ), val( geo ), val( group_nr ), val( group_size ), val( samples ) */
	/*[LG20, [msmc_run:45, spec:uni, geo:pan, group_nr:4, group_size:3, samples:ind1, ind2, ind3], coverage done!, segsites done!]*/
	set val( lg ), val( run ), msmc_gr from msmc_grouping_after_segsites

	output:
	set val( run ), val( lg ), val( msmc_gr.spec ), val( msmc_gr.geo ), val( msmc_gr.group_size ), file( "msmc_run.*.multihetsep.txt" ) into msmc_input_lg

	script:
	"""
	COVDIR="\$BASE_DIR/ressources/coverage_masks/"
	SMP=\$(echo ${msmc_gr.samples}  | \
		sed "s|, |\\n--mask=\${COVDIR}|g; s|^|--mask=\${COVDIR}|g" | \
		sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
		echo \$( cat ) )

	SEGDIR="\$BASE_DIR/2_analysis/msmc/segsites/"
	SEG=\$(echo ${msmc_gr.samples}  | \
		sed "s|, |\\n\${SEGDIR}|g; s|^|\${SEGDIR}|g" | \
		sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
		echo \$( cat ) )

	generate_multihetsep.py \
		\$SMP \
		--mask=\$BASE_DIR/ressources/mappability_masks/${lg}.mapmask.bed.txt.gz \
		--negative_mask=\$BASE_DIR/ressources/indel_masks/indel_mask.${lg}.bed.gz \
		\$SEG > msmc_run.${msmc_gr.run}.${msmc_gr.spec}.${msmc_gr.geo}.${lg}.multihetsep.txt
	"""
}

msmc_input_lg
	.groupTuple()
	.set {msmc_input}

/* run msmc ------------------ */

process msmc_run {
	label "L_190g100h_msmc_run"
	publishDir "../../2_analysis/msmc/output/", mode: 'copy' , pattern: "*.final.txt"
	publishDir "../../2_analysis/msmc/loops/", mode: 'copy' , pattern: "*.loop.txt"

	input:
	set msmc_run, lg , spec, geo, group_size, file( hetsep ) from msmc_input

	output:
	file("*.msmc2.*.txt") into msmc_output

	script:
	"""
	NHAP=\$(echo \$(seq 0 \$((${group_size[0]}*2-1))) | sed 's/ /,/g' )
	INFILES=\$( echo ${hetsep} )

	msmc2 \
		-m 0.00254966 -t 8 \
		-p 1*2+25*1+1*2+1*3 \
		-o run${msmc_run}.${spec[0]}.${geo[0]}.msmc2 \
		-I \${NHAP} \
		\${INFILES}
	"""
}

/* generating MSMC cross coalescence input files (2 inds x 2 species) ----------- */

cc_grouping
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ run:row.run_nr, geo:row.geo, spec_1:row.spec_1, spec_2:row.spec_2, contrast_nr:row.contrast_nr, samples_1:row.samples_1, samples_2:row.samples_2 ] }
	.set { cc_runs }


lg_ch3
	.combine( cc_runs )
	.combine( coverage_cc )
	.combine( segsites_cc )
	.map{[it[0], it[1].run, it[1]]}
	.set{ cc_grouping_after_segsites }

process generate_multihetsep_cc {
	label "L_105g30h_cc_generate_multihetsep"
	publishDir "../../2_analysis/cross_coalescence/input/run_${run}", mode: 'copy' , pattern: "*.multihetsep.txt"
	conda "$HOME/miniconda2/envs/py3"

	input:
	/* content cc_gr: val( run_nr ), val( geo ), val( spec_1 ), val( spec_2 ), val( contrast_nr ), val( samples_1 ), val( samples_2 ) */
	set val( lg ), val( run ), cc_gr from cc_grouping_after_segsites

	output:
	set val( cc_gr.run ), val( lg ), val( cc_gr.spec_1 ), val( cc_gr.spec_2 ), val( cc_gr.geo ), val( cc_gr.contrast_nr ), val( cc_gr.samples_1 ), val( cc_gr.samples_2 ), file( "cc_run.*.multihetsep.txt" ) into cc_input_lg
	/* !! CHECK: hetsep  using ALL samples of species? */
	/* !! CHECK: also - pipe at indel script broken? (filter_indels)  */

	script:
	"""
	COVDIR="\$BASE_DIR/ressources/coverage_masks/"
	SMP1=\$(echo ${cc_gr.samples_1}  | \
		sed "s|, |\\n--mask=\${COVDIR}|g; s|^|--mask=\${COVDIR}|g" | \
		sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
		echo \$( cat ) )
	SMP2=\$(echo ${cc_gr.samples_2}  | \
		sed "s|, |\\n--mask=\${COVDIR}|g; s|^|--mask=\${COVDIR}|g" | \
		sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
		echo \$( cat ) )

	SEGDIR="\$BASE_DIR/2_analysis/msmc/segsites/"
	SEG1=\$(echo ${cc_gr.samples_1}  | \
		sed "s|, |\\n\${SEGDIR}|g; s|^|\${SEGDIR}|g" | \
		sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
		echo \$( cat ) )
	SEG2=\$(echo ${cc_gr.samples_2}  | \
		sed "s|, |\\n\${SEGDIR}|g; s|^|\${SEGDIR}|g" | \
		sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
		echo \$( cat ) )

	generate_multihetsep.py \
		\${SMP1} \
		\${SMP2} \
		--mask=\$BASE_DIR/ressources/mappability_masks/${lg}.mapmask.bed.txt.gz \
		--negative_mask=\$BASE_DIR/ressources/indel_masks/indel_mask.${lg}.bed.gz \
		\${SEG1} \
		\${SEG2} \
		> cc_run.${run}.${cc_gr.spec_1}-${cc_gr.spec_2}.${cc_gr.contrast_nr}.${cc_gr.geo}.${lg}.multihetsep.txt
	"""
}
/*
cc_input_lg
  .groupTuple()
	.set {cc_input}
*/
/* run cross coalescence -------------- */
/*process cc_run {
	label "L_190g30ht24_cc_run"
	publishDir "../../2_analysis/cross_coalescence/output/", mode: 'copy'
	conda "$HOME/miniconda2/envs/py3"

	input:
	set cc_run, lg , spec1, spec2, geo, contr_nr, samples_1, samples_2, file( hetsep ) from cc_input

	output:
	file("cc_run.*.final.txt.gz") into cc_output
*/
	/* !! CHECK: using hetsep index for -I flag? (currently sample ID is used)?
		# --> replaced by index */
/*
	script:
	"""
	INFILES=\$( echo ${hetsep} )
	POP1=\$( echo "${samples_1}" | sed 's/\\[//g; s/, /,/g; s/\\]//g' )
	POP2=\$( echo "${samples_2}" | sed 's/\\[//g; s/, /,/g; s/\\]//g' )

	msmc2 \
		-m 0.00255863 -t 24 \
		-p 1*2+25*1+1*2+1*3 \
		-o cc_run.${cc_run[0]}.${spec1[0]}.msmc \
		-I 0,1,2,3 \ # \${POP1} \
		\${INFILES}

	msmc2 \
		-m 0.00255863 -t 24 \
		-p 1*2+25*1+1*2+1*3 \
		-o cc_run.${cc_run[0]}.${spec2[0]}.msmc \
		-I 4,5,6,7 \ # \${POP2} \
		\${INFILES}

	msmc2 \
		-m 0.00255863 -t 24 \
		-p 1*2+25*1+1*2+1*3 \
		-o cc_run.${cc_run[0]}.cross.msmc \
		-I 0,1,2,3,4,5,6,7 \ # \${POP1},\${POP2} \
		-P 0,0,0,0,1,1,1,1 \
		\${INFILES}

	combineCrossCoal.py \
		cc_run.${cc_run[0]}.cross.msmc.final.txt \
		cc_run.${cc_run[0]}.${spec1[0]}.msmc.final.txt \
		cc_run.${cc_run[0]}.${spec2[0]}.msmc.final.txt | \
		gzip > cc_run.${cc_run[0]}.final.txt.gz
	"""
}
*/