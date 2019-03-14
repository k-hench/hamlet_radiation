#!/usr/bin/env nextflow
/* create channel of linkage groups */
Channel
	.from( 2..12 )
	.into{ admx_ch; admx_loc_ch }

Channel
	.fromFilePairs("1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_locations; vcf_all_samples_pca; vcf_admx; vcf_geno }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

locations_ch
	.combine( vcf_locations )
	.set{ vcf_location_combo }

process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_pca, vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3, vcf_loc_admix )

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

/* 1) PCA section ============== */
/* 1a) PCA (local) -------------- */
process pca_location {
	label "L_20g15h_pca_location"
	publishDir "figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "2_analysis/pca", mode: 'copy' , pattern: "*.gz"
	module "R3.5.2"

	input:
	set val( loc ), file( vcf ), file( pop ) from vcf_loc_pca

	output:
	set file( "${loc}.prime_pca.pdf" ), file( "${loc}.pca.pdf" ), file( "${loc}.exp_var.txt.gz" ), file( "${loc}.scores.txt.gz" ) into pca_loc_out

	script:
	"""
	awk '{print \$1"\\t"\$1}' ${loc}.pop | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > ${loc}.pop.txt

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R ${loc}.pop.txt 6
	"""
}

/* 1b) PCA (global) -------------- */
process pca_all {
	label "L_20g15h_pca_all"
	publishDir "figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "2_analysis/pca", mode: 'copy' , pattern: "*.txt.gz"
	publishDir "1_genotyping/4_phased/", mode: 'copy' , pattern: "*.vcf.gz"
	module "R3.5.2"

	input:
	set vcfId, file( vcf ) from vcf_all_samples_pca

	output:
	set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out
	file( "hamlets_only.vcf.gz*" ) into vcf_hamlets_only

	script:
	"""
	# complete PCA, all samples ------------
	vcfsamplenames ${vcf[0]} | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > all.pop.txt

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R all.pop.txt 6

	# PCA without outgroups ---------------

	vcfsamplenames ${vcf[0]} | \
		grep -v "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" > outgroup.pop

	vcfsamplenames ${vcf[0]} | \
		grep "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > hamlets_only.pop.txt

	vcftools \
		--gzvcf ${vcf[0]} \
		--remove outgroup.pop \
		--recode \
		--stdout | gzip > hamlets_only.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.vcf.gz \$BASE_DIR/R/project_config.R hamlets_only.pop.txt 6

	# PCA without indtgo or gumigutta ---------------

	grep -v "ind\\|gum" hamlets_only.pop.txt > core_hamlets.pop.txt
	cut -f 1 core_hamlets.pop.txt > core_hamlets.pop

	vcftools \
		--gzvcf ${vcf[0]} \
		--keep core_hamlets.pop \
		--recode \
		--stdout | gzip > core_hamlets.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R core_hamlets.vcf.gz \$BASE_DIR/R/project_config.R core_hamlets.pop.txt 6
	"""
}

/* 2) Admixture section ============== */

/* 2a) Admixture (global) -------------- */
process plink12 {
	label 'L_20g2h_plink12'

	input:
	set vcfId, file( vcf ) from vcf_admx

	output:
	set file( "hapmap.ped" ), file( "hapmap.map" ), file( "hapmap.nosex" ), file( "pop.txt" ) into admx_plink
	set file( "GxP_plink.map" ), file( "GxP_plink.ped" ) into plink_GxP

	script:
	"""
	vcfsamplenames ${vcf[0]} | \
		grep -v "tor\\|tab\\|flo" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > pop.txt

	vcftools \
		--gzvcf ${vcf[0]} \
		--plink \
		--out GxP_plink

	plink \
		--file GxP_plink \
		--recode12 \
		--out hapmap
	"""
}

admx_prep  = admx_ch.combine( admx_plink )

process admixture_all {
	label 'L_O88g200h_admixture_all'
	publishDir "2_analysis/admixture/", mode: 'copy' , pattern: "*.Q"

	input:
	set  val( x ), file( ped ), file( map ), file( nosex ), file( pop ) from admx_prep

	output:
	file( "hapmap.all.${x}.P" ) into admx_output
	set val( "dummy" ), file( "log${x}-all.out" ), file( "hapmap.all.${x}.Q" ), file( "pop.${x}.txt" ) into admx_log

	script:
	"""
	mv ${pop} pop.${x}.txt
	admixture --cv ${ped} ${x} | tee log${x}-all.out
	mv hapmap.${x}.P hapmap.all.${x}.P
	mv hapmap.${x}.Q hapmap.all.${x}.Q
	"""
}

process admixture_log {
	label 'L_loc_admixture_log'
	publishDir "2_analysis/admixture/", mode: 'copy' , pattern: "admixture_report*"
	publishDir "figures/admixture", mode: 'copy' , pattern: "*.pdf"
	module "R3.5.2"

	input:
	set val( dummy ), file( logs ), file( admxQ ), file( pop ) from admx_log.groupTuple()

	output:
	file( "admixture_report.txt" ) into admxR_output
	file( "admixture*.pdf" ) into admx_plot_all

	script:
	"""
	grep -h CV log*.out > admixture_report.txt
	Rscript --vanilla \$BASE_DIR/R/plot_admixture.R admixture_report.txt ${pop[0]} \$BASE_DIR/R/project_config.R .all 12
	"""
}

/* 2b) Admixture (local) -------------- */
process plink12_loc {
	label 'L_20g2h_plink12_loc'

	input:
	set val( loc ), file( vcf ), file( pop ) from vcf_loc_admix

	output:
	set val( loc ), file( "hapmap.${loc}.ped" ), file( "hapmap.${loc}.map" ), file( "hapmap.${loc}.nosex" ), file( pop ) into admx_loc_plink

	script:
	"""
	vcftools \
		--gzvcf ${vcf} \
		--plink \
		--out intermediate_plink

	plink \
		--file intermediate_plink \
		--recode12 \
		--out hapmap.${loc}
	"""
}

admx_loc_prep  = admx_loc_ch.combine( admx_loc_plink )

process admixture_loc {
	label 'L_78g10h_admixture_loc'
	publishDir "2_analysis/admixture/${loc}", mode: 'copy' , pattern: "*.Q"

	input:
	set val( x ), val( loc ), file( ped ), file( map ), file( nosex ), file( pop ) from admx_loc_prep

	output:
	set val( loc ), file( "*.P" ) into admxP_loc
	set val( loc ), file( "*.Q" ),file( "log*.out" ), file( "${x}${pop}" ) into admxQL_loc

	script:
	"""
	mv ${pop} ${x}${pop}
	admixture --cv ${ped} ${x} | tee log${x}-${loc}.out
	"""
}

admxQL_loc
	.groupTuple()
	.set { admx_loc_log_sorted }

process admixture_loc_log {
	label 'L_loc_admixture_log_loc'
	publishDir "2_analysis/admixture/${loc}", mode: 'copy' , pattern: "admixture_report*"
	publishDir "figures/admixture", mode: 'copy' , pattern: "*.pdf"
	module "R3.5.2"

	input:
	set val( loc ), file( admxQ ), file( logs ), file( pop ) from admx_loc_log_sorted

	output:
	file( "admixture_report.${loc}.txt" ) into admxR_loc_output
	file( "admixture*.pdf" ) into admx_plot_loc

	script:
	"""
	cat ${pop[0]} | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\2/g' > pop.txt

	grep -h CV log*.out > admixture_report.${loc}.txt
	Rscript --vanilla \$BASE_DIR/R/plot_admixture.R admixture_report.${loc}.txt pop.txt \$BASE_DIR/R/project_config.R .${loc} 8
	"""
	/*input: 1:admixture_reports, 2: sample_ids, 3: project_config, 4:location, 5: max_admixture_K */
}

/* 3) Fst section ============== */
/* Preparation: create all possible species pairs depending on location
   and combine with genotype subset (for the respective location)*/

/* channel content after joinig: set [0:val(loc), 1:file(vcf), 2:file(pop), 3:val(spec1), 4:val(spec2)]*/
bel_pairs_ch = Channel.from( "bel" )
	.join( vcf_loc_pair1 )
	.combine(bel_spec1_ch)
	.combine(bel_spec2_ch)
	.filter{ it[3] < it[5] }
	.map{ it[0,1,2,4,6]}
hon_pairs_ch = Channel.from( "hon" )
	.join( vcf_loc_pair2 )
	.combine(hon_spec1_ch)
	.combine(hon_spec2_ch)
	.filter{ it[3] < it[5] }
	.map{ it[0,1,2,4,6]}
pan_pairs_ch = Channel.from( "pan" )
	.join( vcf_loc_pair3 )
	.combine(pan_spec1_ch)
	.combine(pan_spec2_ch)
	.filter{ it[3] < it[5] }
	.map{ it[0,1,2,4,6]}
bel_pairs_ch.concat( hon_pairs_ch, pan_pairs_ch  ).set { all_fst_pairs_ch }

process fst_run {
	label 'L_32g4h_fst_run'
	publishDir "2_analysis/fst/50k/${loc}", mode: 'copy' , pattern: "*.50k.windowed.weir.fst.gz"
	publishDir "2_analysis/fst/10k/${loc}", mode: 'copy' , pattern: "*.10k.windowed.weir.fst.gz"
	publishDir "2_analysis/fst/logs/${loc}", mode: 'copy' , pattern: "${loc}-${spec1}-${spec2}.log"

	input:
	set val( loc ), file( vcf ), file( pop ), val( spec1 ), val( spec2 ) from all_fst_pairs_ch

	output:
	set val( loc ), file( "*.50k.windowed.weir.fst.gz" ), file( "${loc}-${spec1}-${spec2}.log" ) into fst_50k
	file( "*.10k.windowed.weir.fst.gz" ) into fst_10k_output
	file( "${loc}-${spec1}-${spec2}.log" ) into fst_logs

	script:
	"""
	grep ${spec1} ${pop} > pop1.txt
	grep ${spec2} ${pop} > pop2.txt

	vcftools --gzvcf ${vcf} \
		--weir-fst-pop pop1.txt \
		--weir-fst-pop pop2.txt \
		--fst-window-step 5000 \
		--fst-window-size 50000 \
		--out ${loc}-${spec1}-${spec2}.50k 2> ${loc}-${spec1}-${spec2}.log

	vcftools --gzvcf ${vcf} \
		--weir-fst-pop pop1.txt \
		--weir-fst-pop pop2.txt \
		--fst-window-size 10000 \
		--fst-window-step 1000 \
		--out ${loc}-${spec1}-${spec2}.10k

	gzip *.windowed.weir.fst
	"""
}

/* collect the VCFtools logs to crate a table with the
   genome wide fst values */
process fst_globals {
	label 'L_loc_fst_globals'
	publishDir "2_analysis/fst/logs/", mode: 'copy' , pattern: "fst_globals.txt"
	publishDir "figures/fst", mode: 'copy' , pattern: "global_fst.pdf"
	module "R3.5.2"

	input:
	file( log ) from fst_logs.collect()

	output:
	file( "fst_globals.txt" ) into fst_glob

	script:
	"""
	cat *.log | \
		grep -E 'Weir and Cockerham|--out' | \
		grep -A 3 50k | \
		sed '/^--/d; s/^.*--out //g; s/.50k//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
		paste - - - | \
		cut -f 1,3,5 | \

	sed 's/^\\(...\\)-/\\1\\t/g' > fst_globals.txt
	Rscript --vanilla \$BASE_DIR/R/plot_global_fst.R fst_globals.txt \$BASE_DIR/R/fst_functions.R \$BASE_DIR/R/project_config.R
	"""
}

fst_50k
	.groupTuple()
	.set { fst_50k_sorted }

process plot_fst {
	label 'L_20g2h_plot_fst'
	publishDir "figures/fst", mode: 'copy' , pattern: "*.pdf"
	module "R3.5.2"

	input:
	set val( loc ), file( first_fst ), file( first_log ) from fst_50k_sorted

	output:
	file( "*.pdf" ) into fst_plots

	script:
	"""
	cat *.log | \
		grep -E 'Weir and Cockerham|--out' | \
		grep -A 3 50k | \
		sed '/^--/d; s/^.*--out //g; s/.50k//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
		paste - - - | \
		cut -f 1,3,5 | \
		sed 's/^\\(...\\)-/\\1\\t/g' > fst_${loc}.txt

	Rscript --vanilla \$BASE_DIR/R/plot_fst.R ${loc} fst_${loc}.txt \$BASE_DIR/R/fst_functions.R \$BASE_DIR/R/project_config.R
	"""
}

/* 4) G x P section ============== */
process GxP_run {
	label 'L_20g2h_GxP_binary'

	input:
	set file( map ), file( ped ) from plink_GxP

	output:
	set file( "*.bed" ), file( "*.bim" ),file( "*.fam" ) into plink_binary

	script:
	"""
	# convert genotypes into binary format (bed/bim/fam)
	plink \
		--noweb \
		--file GxP_plink \
		--make-bed \
		--out GxP_plink_binary
	"""
}

Channel
	.fromPath("metadata/phenotypes.sc")
	.set{ phenotypes_raw }

process phenotye_pca {
	label "L_loc_phenotype_pca"
	publishDir "figures/phenotype", mode: 'copy' , pattern: "*.pdf"
	publishDir "2_analysis/phenotype", mode: 'copy' , pattern: "*.gz"
	module "R3.5.2"

	input:
	file( sc ) from phenotypes_raw

	output:
	file( "phenotypes.txt" ) into phenotype_file
	file( "phenotype_pca*.pdf" ) into  phenotype_pca

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/phenotypes_pca.R ${sc} \$BASE_DIR/R/project_config.R
	"""
}

Channel
	.from("Bars", "Lines", "Snout", "Peduncle", "Blue", "Yellow", "Orange", "Tail_transparent","PC1", "PC2", "PC_d1", "abe", "gum", "ind", "may", "nig", "pue", "ran", "uni")
	.set{ traits_ch }

traits_ch.combine( plink_binary ).combine( phenotype_file ).set{ trait_plink_combo }

process gemma_run {
 label 'L_32g4h_GxP_run'
 publishDir "2_analysis/GxP/bySNP/", mode: 'copy'
 module "R3.5.2"

 input:
 set  val( pheno ), file( bed ), file( bim ), file( fam ), file( pheno_file ) from trait_plink_combo

 output:
 file("*.GxP.txt.gz") into gemma_results

 script:
	"""
	source \$BASE_DIR/sh/body.sh
	BASE_NAME=\$(echo  ${fam} | sed 's/.fam//g')

	mv ${fam} \$BASE_NAME-old.fam
	cp \${BASE_NAME}-old.fam ${fam}

	# 1) replace the phenotype values
	Rscript --vanilla \$BASE_DIR/R/assign_phenotypes.R ${fam} ${pheno_file} ${pheno}

	# 2) create relatedness matrix of samples using gemma
	gemma -bfile \$BASE_NAME -gk 1 -o ${pheno}

	# 3) fit linear model using gemma (-lm)
	gemma -bfile \$BASE_NAME -lm 4 -miss 0.1 -notsnp -o ${pheno}.lm

	# 4) fit linear mixed model using gemma (-lmm)
	gemma -bfile \$BASE_NAME -k output/${pheno}.cXX.txt -lmm 4 -o ${pheno}.lmm

	# 5) reformat output
	sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/${pheno}.lm.assoc.txt | \
		cut -f 2,3,9-14 | body sort -k1,1 -k2,2n | gzip > ${pheno}.lm.GxP.txt.gz
	sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/${pheno}.lmm.assoc.txt | \
		cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > ${pheno}.lmm.GxP.txt.gz
	"""
}

Channel
	.from([[50000, 5000], [10000, 1000]])
	.set{ gxp_smoothing_levels }

gemma_results.combine( gxp_smoothing_levels ).set{ gxp_smoothing_input }

process gemma_smooth {
	label 'L_20g2h_GxP_smooth'
	publishDir "2_analysis/GxP/", mode: 'copy'

	input:
	set file( lm ), file( lmm ), val( win ), val( step ) from gxp_smoothing_input

	output:
	set val( win ), file( "*.lm.*k.txt.gz" ) into gxp_lm_smoothing_output
	set val( win ), file( "*.lmm.*k.txt.gz" ) into gxp_lmm_smoothing_output

	script:
	"""
	\$BASE_DIR/sh/gxp_slider ${lm} ${win} ${step}
	\$BASE_DIR/sh/gxp_slider ${lmm} ${win} ${step}
	"""
}

gxp_lm_smoothing_output
	.groupTuple()
	.filter{ it[0] == 50000 }
	.subscribe{ println it }
/*	.set{ grouped_smoothing_output }*/

/*
process gemma_plot {
	label 'L_loc_GxP_plot'
	publishDir "figures/gxp", mode: 'copy' , pattern: "*.pdf"
	module "R3.5.2"

	input:
	set file( win ), file( gxp ) from grouped_smoothing_output

	output:
	file( "*.pdf" ) into gxp_plot_output

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/plot_gxp.R \$BASE_DIR/R/project_config.R
	"""
} */