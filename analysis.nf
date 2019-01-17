#!/usr/bin/env nextflow
Channel
	.from( 2..15 )
	.into{ admx_ch; admx_loc_ch }

Channel
	.fromFilePairs("1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_phylo; vcf_locations; vcf_all_samples_pca; vcf_admx; vcf_geno }

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
			vcfsamplenames phased_mac2.vcf.gz | \
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

 script:
 """
 vcfsamplenames ${vcf[0]} | \
 		grep -v "tor\\|tab\\|flo" | \
		awk '{print \$1"\\t"\$1}' | \
		sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > pop.txt


 vcftools \
 		--gzvcf ${vcf[0]} \
		--plink \
		--out intermediate_plink

 plink \
 		--file intermediate_plink \
		--recode12 \
		--out hapmap
 """
}

admx_prep  = admx_ch.combine( admx_plink )

process admixture_all {
    label 'L_88g30h_admixture_all'
    publishDir "2_analysis/admixture/", mode: 'copy' , pattern: "*.Q"

    input:
    set val( x ), file( ped ), file( map ), file( nosex ), file( pop ) from admx_prep

    output:
    file( "hapmap.all.${x}.P" ) into admx_output
		set file( "log${x}-all.out" ), file( "hapmap.all.${x}.Q" ), file( pop ) into admx_log

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
  set file( logs ), file( admxQ ), file( pop ) from admx_log.collect()

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
  label 'L_190g4h_fasttree'

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
  label 'L_300g30h_fasttree'
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
	"""cat *.log | \
  grep -E 'Weir and Cockerham|--out' | \
  grep -A 3 50k | \
  sed '/^--/d; s/^.*--out //g; s/.50k//g; /^Output/d; s/Weir and Cockerham //g; s/ Fst estimate: /\t/g' | \
  paste - - - | \
  cut -f 1,3,5 | \
  sed 's/^\\(...\\)-/\\1\\t/g' > fst_${loc}.txt

	Rscript --vanilla \$BASE_DIR/R/plot_fst.R ${loc} fst_${loc}.txt \$BASE_DIR/R/fst_functions.R \$BASE_DIR/R/project_config.R
	"""
}

/* 5) Twisst section ============== */
process vcf2geno_loc {
  label 'L_20g15h_vcf2geno'

  input:
	set val( loc ), file( vcf ), file( pop ) from vcf_loc_twisst

  output:
  set val(loc), file( "${loc}.geno.gz" ), file( pop ) into snp_geno_twisst

  script:
  """
  python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
    -i ${vcf} | gzip > ${loc}.geno.gz
  """
}

Channel.from( 50, 100, 200 ).set{ twisst_window_types }
snp_geno_twisst.combine( twisst_window_types ).set{ twisst_input_ch }

process twisst_prep {
  label 'L_120g30h6t_prep_twisst'

  input:
  set val( loc ), file( geno ), file( pop ), val( twisst_w ) from twisst_input_ch

	output:
	set val( loc ), file( geno ), file( pop ), val( twisst_w ), file( "${loc}.w${twisst_w}.phyml_bionj.trees.gz" ) into twisst_prep_ch

  script:
   """
   python \$SFTWR/genomics_general/phylo/phyml_sliding_windows.py \
      -g ${geno} \
      --windType sites \
      -w ${twisst_w} \
      --prefix ${loc}.w${twisst_w}.phyml_bionj \
      --model HKY85 \
      --optimise n \
      --threads 6
	 """
}

process twisst_run {
  label 'L_120g30h6t_run_twisst'
  publishDir "2_analysis/twisst/", mode: 'copy'

  input:
	set val( loc ), file( geno ), file( pop ), val( twisst_w ), file( tree ) from twisst_prep_ch

	output:
  set val( loc ), val( twisst_w ), file( "*.weights.tsv.gz" ), file( "*.data.tsv" ) into ( twisst_output )

  script:
   """
	awk '{print \$1"\\t"\$1}' ${pop} | \
		sed 's/\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' | \
		cut -f 1,3 | \
		awk '{print \$1"_A\\t"\$2"\\n"\$1"_B\\t"\$2}' > ${loc}.twisst_pop.txt

	TWISST_POPS=\$( cut -f 2 ${loc}.twisst_pop.txt | sort | uniq | echo \$( cat ) | sed 's/ / -g /g; s/^/-g /' )

   python \$SFTWR/twisst/run_twisst_parallel.py \
      --method complete \
      --threads 6 \
      -t ${tree} \
      \$TWISST_POPS \
      --groupsFile ${loc}.twisst_pop.txt | \
      gzip > ${loc}.w${twisst_w}.phyml_bionj.weights.tsv.gz
    """
  }