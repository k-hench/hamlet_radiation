#!/usr/bin/env nextflow
/* create channel of linkage groups */
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.map{ "LG" + it }
	.into{ lg_fasttree; lg_twisst }

Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations;  vcf_geno; vcf_fasttree_whg }

Channel
	.from( "bel", "hon" )
	.set{ locations_ch }

locations_ch
	.combine( vcf_locations )
	.set{ vcf_location_combo }

Channel
	.from( "all", "bel", "hon", "pan" )
	.set{ locations4_ch }

Channel
	.from( "whg_no_og", "no_musks" )
	.set{ whg_modes }

locations4_ch
	.combine( vcf_fasttree_whg )
	.combine( whg_modes )
	.set{ vcf_fasttree_whg_location_combo }

process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_twisst )

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
		--stdout | gzip > ${loc}.vcf.gz
	"""
}


/* 1) fasttree section ============== */
/* 1.1) --- by LG --- */
vcf_geno.combine( lg_fasttree ).into{ fasttree_geno; fasttree_no_og_geno  }

process vcf2geno_all {
	label 'L_20g15h_vcf2geno_all'

	input:
	set vcfId, file( vcf ), val( lg ) from fasttree_geno

	output:
	set val( "all" ), val( lg ), file( "output.all.${lg}.geno.gz" ) into snp_geno_tree_all

	script:
	"""
	vcftools \
	--gzvcf ${vcf[0]} \
	--chr ${lg} \
	--recode \
	--stdout | gzip > intermediate.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i  intermediate.vcf.gz | gzip > output.all.${lg}.geno.gz

	rm intermediate.vcf.gz
	"""
}

process vcf2geno_no_og {
	label 'L_20g15h_vcf2geno_no_og'

	input:
	set vcfId, file( vcf ), val( lg ) from fasttree_no_og_geno

	output:
	set val( "no_og" ), val( lg ), file( "output.no_og.${lg}.geno.gz" ) into snp_geno_tree_no_og

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep "tor\\|tab" > og.pop

	vcftools \
	--gzvcf ${vcf[0]} \
	--chr ${lg} \
	--remove og.pop \
	--recode \
	--stdout | gzip > intermediate.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i  intermediate.vcf.gz | gzip > output.no_og.${lg}.geno.gz

	rm intermediate.vcf.gz
	"""
}

snp_geno_tree_all
	.concat( snp_geno_tree_no_og )
	.set{ snp_geno_tree }

process fasttree_prep {
	label 'L_190g15h_fasttree_prep'

	input:
	set val( type ), val( lg ), file( geno ) from snp_geno_tree

	output:
	set val( type ), val( lg ), file( "all_samples.${type}.${lg}.SNP.fa" ) into ( fasttree_prep_ch )

	script:
	"""
	python \$SFTWR/genomics_general/genoToSeq.py -g ${geno} \
		-s  all_samples.${type}.${lg}.SNP.fa \
		-f fasta \
		--splitPhased
	"""
}

process fasttree_run {
	label 'L_190g100h_fasttree_run'
	publishDir "../../2_analysis/fasttree/", mode: 'copy'

	input:
	set val( type ), val( lg ), file( fa ) from fasttree_prep_ch

	output:
	file( "all_samples.${type}.${lg}.SNP.tree" ) into ( fasttree_output )

	script:
	"""
	fasttree -nt ${fa} > all_samples.${type}.${lg}.SNP.tree
	"""
}

/* 1.2) --- whole genome --- */

process subset_vcf_by_location_whg {
	label "L_20g2h_subset_vcf_whg"

	input:
	set val( loc ), vcfId, file( vcf ), val( mode ) from vcf_fasttree_whg_location_combo

	output:
	set val( mode ), val( loc ), file( "${loc}.${mode}.whg.geno.gz" ) into snp_geno_tree_whg_no_og

	script:
	"""
	if [ "${loc}" == "all" ];then
		vcfsamplenames ${vcf[0]} | \
			grep -v tor | \
			grep -v tab > ${loc}.pop
	else
		vcfsamplenames ${vcf[0]} | \
			grep ${loc} | \
			grep -v tor | \
			grep -v tab > ${loc}.pop
	fi

	if [ "${mode}" == "no_musks" ];then
		DROP_CHRS="--not-chr LG04 --not-chr LG07 --not-chr LG08 --not-chr LG09 --not-chr LG12 --not-chr LG17 --not-chr LG23"
	fi

	vcftools --gzvcf ${vcf[0]} \
		--keep ${loc}.pop \
		\$DROP_CHRS \
		--mac 3 \
		--recode \
		--stdout | gzip > ${loc}.${mode}.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i  ${loc}.${mode}.vcf.gz | gzip > ${loc}.${mode}.whg.geno.gz
	"""
}

process fasttree_whg_prep {
	label 'L_190g4h_fasttree_whg_prep'

	input:
	set val( mode ), val( loc ), file( geno ) from snp_geno_tree_whg_no_og

	output:
	set val( mode ), val( loc ), file( "all_samples.${loc}.${mode}.whg.SNP.fa" ) into ( fasttree_whg_prep_ch )

	script:
	"""
	python \$SFTWR/genomics_general/genoToSeq.py -g ${geno} \
		-s  all_samples.${loc}.${mode}.whg.SNP.fa \
		-f fasta \
		--splitPhased
	"""
}

process fasttree_whg_run {
	label 'L_190g100h_fasttree_run'
	publishDir "../../2_analysis/fasttree/", mode: 'copy'

	input:
	set val( mode ), val( loc ), file( fa ) from fasttree_whg_prep_ch

	output:
	file( "all_samples.${loc}.${mode}.whg.SNP.tree" ) into ( fasttree_whg_output )

	script:
	"""
	fasttree -nt ${fa} > all_samples.${loc}.${mode}.whg.SNP.tree
	"""
}

Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd_bi-allelic.mito.vcf.{gz,gz.tbi}")
	.set{ vcf_mito }

process fasttree_mito {
	label 'L_20g2h_fasttree_mito'
	publishDir "../../2_analysis/fasttree/", mode: 'copy'

	input:
	set val( vcf_id ), file( vcf ) from vcf_mito

	output:
	file( "only_a.mito.SNP.tree" ) into ( fasttree_mito_output )

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
	grep -v tor | \
	grep -v tab > vcf.pop

	vcftools --gzvcf ${vcf[0]} \
		--keep vcf.pop \
		--mac 2 \
		--recode \
		--stdout | gzip > mito.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i  mito.vcf.gz | \
		gzip > mito.geno.gz

	python \$SFTWR/genomics_general/genoToSeq.py \
		-g mito.geno.gz \
		-s  all_samples.mito.fa \
		-f fasta \
		--splitPhased

	samtools faidx all_samples.mito.fa

	grep "_A" all_samples.mito.fa.fai | \
		awk -v OFS='\t' '{print \$1,"0",\$2,\$1}' | \
		fastaFromBed -fi all_samples.mito.fa -bed stdin -name -fo only_a.mito.fa

	fasttree -nt only_a.mito.fa > only_a.mito.SNP.tree
	"""
}
/*--------- tree construction -----------*/
/*
process plot_tree {
	label '32g1h.fasttree_plot'
	publishDir "../../out/fasttree/", mode: 'symlink'
	module "R3.5.2"

	input:
	file( tree ) from fasttree_output

	output:
	file( "*.pdf" ) into fasttree_plot

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/plot_tree.R ${tree} \$BASE_DIR/vcf_samples.txt
	"""
}
*/

/* 2) Twisst section ============== */

/* MUTE:
vcf_loc_twisst
	.combine( lg_twisst )
	.set{ vcf_loc_lg_twisst }
*/

/* MUTE: python thread conflict - run locally and feed into ressources/plugin
process vcf2geno_loc {
	label 'L_20g15h_vcf2geno'

	input:
	set val( loc ), file( vcf ), file( pop ), val( lg ) from vcf_loc_lg_twisst

	output:
	set val( loc ), val( lg ), file( "${loc}.${lg}.geno.gz" ), file( pop ) into snp_geno_twisst

	script:
	"""
	vcftools \
		--gzvcf ${vcf} \
		--chr ${lg} \
		--recode \
		--stdout |
		gzip > intermediate.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
	  -i intermediate.vcf.gz | gzip > ${loc}.${lg}.geno.gz
	"""
}
*/
/* MUTE: python thread conflict - run locally and feed into ressources/plugin
Channel.from( 50 ).set{ twisst_window_types }

snp_geno_twisst.combine( twisst_window_types ).set{ twisst_input_ch }
*/
/*
process twisst_prep {
  label 'L_G120g40h_prep_twisst'

  input:
  set val( loc ), val( lg ), file( geno ), file( pop ), val( twisst_w ) from twisst_input_ch.filter { it[0] != 'pan' }

	output:
	set val( loc ), val( lg ), file( geno ), file( pop ), val( twisst_w ), file( "*.trees.gz" ), file( "*.data.tsv" ) into twisst_prep_ch

  script:
   """
	module load intel17.0.4 intelmpi17.0.4

   mpirun \$NQSII_MPIOPTS -np 1 \
	python \$SFTWR/genomics_general/phylo/phyml_sliding_windows.py \
      -g ${geno} \
      --windType sites \
      -w ${twisst_w} \
      --prefix ${loc}.${lg}.w${twisst_w}.phyml_bionj \
      --model HKY85 \
      --optimise n \
		--threads 1
	 """
}
*/
/* MUTE: python thread conflict - run locally and feed into ressources/plugin
process twisst_run {
	label 'L_G120g40h_run_twisst'
	publishDir "../../2_analysis/twisst/", mode: 'copy'

	input:
	set val( loc ), val( lg ), file( geno ), file( pop ), val( twisst_w ), file( tree ), file( data ) from twisst_prep_ch

	output:
	set val( loc ), val( lg ), val( twisst_w ), file( "*.weights.tsv.gz" ), file( "*.data.tsv" ) into ( twisst_output )

	script:
	"""
	module load intel17.0.4 intelmpi17.0.4

	awk '{print \$1"\\t"\$1}' ${pop} | \
	sed 's/\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' | \
	cut -f 1,3 | \
	awk '{print \$1"_A\\t"\$2"\\n"\$1"_B\\t"\$2}' > ${loc}.${lg}.twisst_pop.txt

	TWISST_POPS=\$( cut -f 2 ${loc}.${lg}.twisst_pop.txt | sort | uniq | paste -s -d',' | sed 's/,/ -g /g; s/^/-g /' )

	mpirun \$NQSII_MPIOPTS -np 1 \
	python \$SFTWR/twisst/twisst.py \
	  --method complete \
	  -t ${tree} \
	  -T 1 \
	  \$TWISST_POPS \
	  --groupsFile ${loc}.${lg}.twisst_pop.txt | \
	  gzip > ${loc}.${lg}.w${twisst_w}.phyml_bionj.weights.tsv.gz
	"""
}
*/

process twisst_plugin {
	label 'L_G120g40h_twisst_plugin'
	publishDir "../../2_analysis/twisst/", mode: 'copy'

	input:
	set val( loc ), file( vcf ), file( pop ), val( lg ) from vcf_loc_twisst.combine( lg_twisst )

	output:
	set val( loc ), val( lg ), file( "*.weights.tsv.gz" ) into ( twisst_output )

	script:
	"""
	module load intel17.0.4 intelmpi17.0.4

	awk '{print \$1"\\t"\$1}' ${pop} | \
	sed 's/\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' | \
	cut -f 1,3 | \
	awk '{print \$1"_A\\t"\$2"\\n"\$1"_B\\t"\$2}' > ${loc}.${lg}.twisst_pop.txt

	TWISST_POPS=\$( cut -f 2 ${loc}.${lg}.twisst_pop.txt | sort | uniq | paste -s -d',' | sed 's/,/ -g /g; s/^/-g /' )

	mpirun \$NQSII_MPIOPTS -np 1 \
	python \$SFTWR/twisst/twisst.py \
	  --method complete \
	  -t \$BASE_DIR/ressources/plugin/trees/${loc}/${loc}.${lg}.w50.phyml_bionj.trees.gz \
	  \$TWISST_POPS \
	  --groupsFile ${loc}.${lg}.twisst_pop.txt | \
	  gzip > ${loc}.${lg}.w50.phyml_bionj.weights.tsv.gz
	"""
}