#!/usr/bin/env nextflow
/* create channel of linkage groups */
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.map{ "LG" + it }
	.into{ lg_ch1; lg_ch2; lg_ch3; lg_twisst }

Channel
	.from( 2..15 )
	.into{ admx_ch; admx_loc_ch }

Channel
	.fromFilePairs("1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_phylo; vcf_locations; vcf_all_samples_pca; vcf_admx; vcf_geno; vcf_msmc }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

locations_ch.combine( vcf_locations ).set{ vcf_location_combo }

process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
	set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	output:
	set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_pca, vcf_loc_pair1, vcf_loc_pair2, vcf_loc_pair3, vcf_loc_admix, vcf_loc_twisst  )

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
		--stdout | bgzip > ${loc}.vcf.gz
	"""
}

/* 1) PCA section ============== */
/* 1a) PCA (local) -------------- */
process pca_location {
	label "L_20g15h_pca_location"
	publishDir "figures/pca", mode: 'copy' , pattern: "*.pdf"
	publishDir "2_analysis/pca", mode: 'copy' , pattern: "*.gz"

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
		--stdout | bgzip > hamlets_only.vcf.gz

	tabix hamlets_only.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.vcf.gz \$BASE_DIR/R/project_config.R hamlets_only.pop.txt 6

	# PCA without indtgo or gumigutta ---------------

	grep -v "ind\\|gum" hamlets_only.pop.txt > oscar_special.pop.txt
	cut -f 1 oscar_special.pop.txt > oscar_special.pop

	vcftools \
		--gzvcf ${vcf[0]} \
		--keep oscar_special.pop \
		--recode \
		--stdout | bgzip > oscar_special.vcf.gz

	tabix oscar_special.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R oscar_special.vcf.gz \$BASE_DIR/R/project_config.R oscar_special.pop.txt 6
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
	set vall( "dummy" ), file( "log${x}-all.out" ), file( "hapmap.all.${x}.Q" ), file( pop ) into admx_log

	script:
	"""
	admixture --cv ${ped} ${x} | tee log${x}-all.out
	mv hapmap.${x}.P hapmap.all.${x}.P
	mv hapmap.${x}.Q hapmap.all.${x}.Q
	"""
}

process admixture_log {
	label 'L_loc_admixture_log'
	publishDir "2_analysis/admixture/", mode: 'copy' , pattern: "admixture_report*"
	publishDir "figures/admixture", mode: 'copy' , pattern: "*.pdf"

	input:
	set vall( dummy ), file( logs ), file( admxQ ), file( pop ) from admx_log.groupTuple()

	output:
	file( "admixture_report.txt" ) into admxR_output
	file( "admixture*.pdf" ) into admx_plot_all

	script:
	"""
	grep -h CV log*.out > admixture_report.txt
	Rscript --vanilla \$BASE_DIR/R/plot_admixture.R admixture_report.txt pop.txt \$BASE_DIR/R/project_config.R .all 8
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

/* 3) fasttree section ============== */
process vcf2geno {
	label 'L_20g15h_vcf2geno'

	input:
	set vcfId, file( vcf ) from vcf_geno

	output:
	file( "output.geno.gz" ) into snp_geno_tree /*geno_output*/

	script:
	"""
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i ${vcf[0]} | gzip > output.geno.gz
	"""
}
/* Section removed (script doesn't finish, alsoSNPs
   filtering condition should allready be fullfilled)*/
/*
process geno_snp {
	label 'L_32g4h4t_geno_snp'

	input:
	file( geno ) from geno_output

	output:
	file( "output.geno.SNP.gz" ) into ( snp_geno_twisst, snp_gene_tree )

	script:
	"""
	python \$SFTWR/genomics_general/filterGenotypes.py \
		-i ${geno} \
		--minAlleles 2 \
		-o output.geno.SNP.gz \
		--threads 4
 """
}*/

process fasttree_prep {
	label 'L_190g15h_fasttree_prep'

	input:
	file( geno ) from snp_geno_tree

	output:
	file( "all_samples.SNP.fa" ) into ( fasttree_prep_ch )

	script:
	"""
	python \$SFTWR/genomics_general/genoToSeq.py -g ${geno} \
		-s  all_samples.SNP.fa \
		-f fasta \
		--splitPhased
	"""
}

process fasttree_run {
	label 'L_300g99h_fasttree_run'
	publishDir "2_analysis/fasttree/", mode: 'copy'

	input:
	file( fa ) from fasttree_prep_ch

	output:
	file( " all_samples.SNP.tree" ) into ( fasttree_output )

	script:
	"""
	fasttree -nt ${fa} > all_samples.SNP.tree
	"""
}
/*--------- tree construction -----------*/
/*
process plot_tree {
	label '32g1h.fasttree_plot'
	publishDir "out/fasttree/", mode: 'symlink'

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

/* 4) Fst section ============== */
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
	publishDir "figures/fst", mode: 'copy' , pattern: "*.png"

	input:
	set val( loc ), file( first_fst ), file( first_log ) from fst_50k_sorted

	output:
	file( "*.png" ) into fst_plots

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

vcf_loc_twisst
	.combine( lg_twisst )
	.set{ vcf_loc_lg_twisst }


/* 5) Twisst section ============== */
process vcf2geno_loc {
	label 'L_20g15h_vcf2geno'

	input:
	set val( loc ), file( vcf ), file( pop ), val( lg )  from vcf_loc_lg_twisst

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

Channel.from( 50 ).set{ twisst_window_types }

snp_geno_twisst.combine( twisst_window_types ).set{ twisst_input_ch }
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
/*
process twisst_run {
	label 'L_G120g40h_run_twisst'
	publishDir "2_analysis/twisst/", mode: 'copy'

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
/* 6) G x P section ============== */
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
	file( "*k.txt.gz" ) into gxp_smoothing_output

	script:
	"""
	\$BASE_DIR/sh/gxp_slider ${lm} ${win} ${step}
	\$BASE_DIR/sh/gxp_slider ${lmm} ${win} ${step}
	"""
}

/* 7) Msmc section ============== */

Channel
	.fromFilePairs("1_genotyping/3_gatk_filtered/filterd_bi-allelic.vcf.{gz,gz.tbi}")
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
	.set { depth_by_sample_ch }

/* create channel from bam files and add sample id */
Channel
	.fromPath( '1_genotyping/0_dedup_bams/*.bam' )
	.map{ file ->
				def key = file.name.toString().tokenize('.').get(0)
				return tuple(key, file)}
				.set{ sample_bams }

/* combine sample bams and sequencing depth */
sample_bams
	.join(depth_by_sample_ch)
	.set{ sample_bam_and_depth }

/* multiply the sample channel by the linkage groups */
sample_bam_and_depth
	.combine( vcf_msmc )
	.combine( lg_ch1 )
	.set{ samples_msmc }

/* split vcf by individual ----------------------------- */
process split_vcf_by_individual {
	label 'L_20g15h_split_by_vcf'

	input:
	set val( id ), file( bam ), val( sites ), val( depth ), file( vcf ), val( lg ) from samples_msmc

	output:
	set val( id ), val( lg ), file( bam ), val( depth ), file( "phased_mac2.${id}.${L}.vcf.gz" ) into ( sample_vcf, sample_vcf2 )

	script:
	"""
	gatk --java-options "-Xmx10G"
		SelectVariants \
		-R \$REF_GENOME \
		-V ${vcf} \
		-sn ${id} \
		-L ${lg}\
		-o phased_mac2.${id}.${L}.vcf.gz
	"""
}

process bam_caller {
	label 'L_36g47h_bam_caller'
	publishDir "ressources/coverage_masks", mode: 'copy' , pattern: "*.coverage_mask.bed.gz"
	conda "$HOME/miniconda2/envs/py3"

	input:
	set val( id ), val( lg ), file( bam ), val( depth ), file( vcf ) from sample_vcf

	output:
	set val( id ), val( lg ), file( "*.bam_caller.vcf.gz" ), file( "*.coverage_mask.bed.gz" ) into coverage_by_sample_lg

	script:
	"""
	samtools mpileup -q 25 -Q 20 -C 50 -u -r ${lg} -f \$REF_GENOME ${bam} | \
		bcftools call -c -V indels | \
		\$BASE_DIR/py/bamHamletCaller.py ${depth} ${ID}.${LG}.coverage_mask.bed.gz | \
		gzip -c > ${ID}.${LG}.bam_caller.vcf.gz
	"""
}

process generate_segsites {
	label "L_36g47h_msmc_generate_segsites"
	publishDir "2_analysis/msmc/segsites", mode: 'copy' , pattern: "*.covered_sites.bed.txt.gz"
	conda "$HOME/miniconda2/envs/py3"

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
	publishDir "2_analysis/msmc/setup", mode: 'copy'

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
	.map{ row -> [ msmc_run:row.msmc_run, spec:row.spec, geo:row.geo, group_nr:row.group_nr, group_size:row.group_size, samples:row.samples ] }
	.set { msmc_runs }

/* wait for bam_caller and generate_segsites to finish: */

coverage_by_sample_lg.collect().map{ [ it ] }.into{ coverage_done; coverage_cc }
segsites_by_sample_lg.collect().map{ [ it ] }.into{ segsites_done; segsites_cc }

lg_ch2
	.combine( msmc_runs )
	.combine( coverage_done )
	.combine( segsites_done )
	.set{ msmc_grouping_after_segsites }

/* generating MSMC input files (4 inds per species) ----------- */
/*
process generate_multihetsep {
	label "L_120g40h_msmc_generate_multihetsep"
	publishDir "2_analysis/msmc/input/run_${msmc_gr.msmc_run}", mode: 'copy' , pattern "*.multihetsep.txt"
	conda "$HOME/miniconda2/envs/py3"

	input:
	/* content msmc_gr: val( msmc_run ), val( spec ), val( geo ), val( group_nr ), val( group_size ), val( samples ) */
	/*
	set val( lg ), msmc_gr, coverage, segsites from msmc_grouping_after_segsites

	output:
	set val( msmc_gr.msmc_run ), val( lg ), val( msmc_gr.spec ), val( msmc_gr.geo ), val( msmc_gr.group_size ), file( "msmc_run.*.multihetsep.txt" ) into msmc_input_lg

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
		sed "s/\$/.${lg}.covered_sites.bed.txt.gz/g" | \
		echo \$( cat ) )

	generate_multihetsep.py \
		\$SMP \ # (--mask=\$COV/{sample}.{lg}.coverage_mask.bed.gz ...)
		--mask=\$BASE_DIR/ressources/mappability_masks/${lg}.mapmask.bed.txt.gz \
		--negative_mask=\$BASE_DIR/ressources/indel_masks/indel_mask.${lg}.bed.gz \
		\$SEG \ # (\${SEGDIR}/{sample}.{lg}.covered_sites.bed.txt.gz ...)
		> msmc_run.${msmc_gr.msmc_run}.${msmc_gr.spec}.${msmc_gr.geo}.${lg}.multihetsep.txt
	"""
}

msmc_input_lg
	.groupTuple()
	.set {msmc_input}
*/
/* run msmc ------------------ */
/*
process msmc_run {
	label "L_190g100h_msmc_run"
	publishDir "2_analysis/msmc/output/", mode: 'copy'

	input:
	set msmc_run, lg , spec, geo, group_size, file( hetsep ) from msmc_input

	output:
	file("*.msmc2") into msmc_output

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
*/
/* generating MSMC cross coalescence input files (2 inds x 2 species) ----------- */
/*
cc_grouping
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ run_nr:row.run_nr, geo:row.geo, spec_1:row.spec_1, spec_2:row.spec_2, contrast_nr:row.contrast_nr, samples_1:row.samples_1, samples_2:row.samples_2 ] }
	.set { cc_runs }

lg_ch3
	.combine( cc_runs )
	.combine( coverage_cc )
	.combine( segsites_cc )
	.set{ cc_grouping_after_segsites }
*/
/*
process generate_multihetsep_cc {
	label "L_105g30h_cc_generate_multihetsep"
	publishDir "2_analysis/cross_coalescence/input/run_${cc_gr.run_nr}", mode: 'copy' , pattern "*.multihetsep.txt"
	conda "$HOME/miniconda2/envs/py3"

	input:
*/
	/* content cc_gr: val( run_nr ), val( geo ), val( spec_1 ), val( spec_2 ), val( contrast_nr ), val( samples_1 ), val( samples_2 ) */
/*
	set val( lg ), cc_gr, coverage, segsites from cc_grouping_after_segsites

	output:
	set val( cc_gr.run_nr ), val( lg ), val( cc_gr.spec_1 ), val( cc_gr.spec_2 ), val( cc_gr.geo ), val( cc_gr.contrast_nr ), val( cc_gr.samples_1 ), val( cc_gr.samples_2 ), file( "cc_run.*.multihetsep.txt" ) into cc_input_lg
*/
	/* !! CHECK: hetsep  using ALL samples of species? */
	/* !! CHECK: also - pipe at indel script broken? (filter_indels)  */
/*
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
		sed "s/\$/.${lg}.covered_sites.bed.txt.gz/g" | \
		echo \$( cat ) )
	SEG2=\$(echo ${cc_gr.samples_2}  | \
		sed "s|, |\\n\${SEGDIR}|g; s|^|\${SEGDIR}|g" | \
		sed "s/\$/.${lg}.covered_sites.bed.txt.gz/g" | \
		echo \$( cat ) )

	generate_multihetsep.py \
		\${SMP1} \
		\${SMP2} \
		--mask=\$BASE_DIR/ressources/mappability_masks/${lg}.mapmask.bed.txt.gz \
		--negative_mask=\$BASE_DIR/ressources/indel_masks/indel_mask.${lg}.bed.gz \
		\${SEG1} \
		\${SEG2} \
		> cc_run.${cc_gr.run_nr}.${cc_gr.spec_1}-${cc_gr.spec_2}.${cc_gr.contrast_nr}.${cc_gr.geo}.${lg}.multihetsep.txt
	"""
}

cc_input_lg
  .groupTuple()
	.set {cc_input}
*/
/* run cross coalescence -------------- */
/*process cc_run {
	label "L_190g30ht24_cc_run"
	publishDir "2_analysis/cross_coalescence/output/", mode: 'copy'
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