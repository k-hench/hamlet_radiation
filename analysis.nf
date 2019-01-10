#!/usr/bin/env nextflow
Channel
	.from( 1..15 )
	.into{ admx_ch; admx_loc_ch }

Channel
	.fromFilePairs("1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_phylo; vcf_locations; vcf_all_samples_pca; vcf_admx; vcf_geno }

Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

vcf_location_combo = locations_ch.combine( vcf_locations )

process subset_vcf_by_location {
	   label "L_20g2h_subset_vcf"

	   input:
		 set val( loc ), vcfId, file( vcf ) from vcf_location_combo

	   output:
	   set val( loc ), file( "${loc}.vcf.gz" ), file( "${loc}.pop" ) into ( vcf_loc_pca, vcf_loc_fst, vcf_loc_admix )

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
		publishDir "figures/pca", mode: 'move' , pattern: "*.pdf"
		publishDir "2_analysis/pca", mode: 'move' , pattern: "*.gz"

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
		publishDir "figures/pca", mode: 'move' , pattern: "*.pdf"
		publishDir "2_analysis/pca", mode: 'move' , pattern: "*.txt.gz"
		publishDir "1_genotyping/4_phased/", mode: 'move' , pattern: "*.vcf.gz"

		input:
		set vcfId, file( vcf ) from vcf_all_samples_pca

		output:
		set file( "*.prime_pca.pdf" ), file( "*.pca.pdf" ), file( "*.exp_var.txt.gz" ), file( "*.scores.txt.gz" ) into pca_all_out
		file( "hamlets_only.vcf.gz*" ) into vcf_hamlets_only

		script:
		"""
		vcfsamplenames ${vcf[0]} | \
			awk '{print \$1"\\t"\$1}' | \
			sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > all.pop.txt

		Rscript --vanilla \$BASE_DIR/R/vcf2pca.R ${vcf[0]} \$BASE_DIR/R/project_config.R all.pop.txt 6

		vcfsamplenames phased_mac2.vcf.gz | \
			grep -v "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" > outgroup.pop

		vcfsamplenames phased_mac2.vcf.gz | \
			grep -v "abe\\|gum\\|ind\\|may\\|nig\\|pue\\|ran\\|uni" | \
				awk '{print \$1"\\t"\$1}' | \
				sed 's/\\t.*\\(...\\)\\(...\\)\$/\\t\\1\\t\\2/g' > hamlets_only.pop.txt

		vcftools \
			--gzvcf ${vcf[0]} \
			--remove outgroup.pop \
			--recode \
			--stdout | bgzip > hamlets_only.vcf.gz

		tabix hamlets_only.vcf.gz

		Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.vcf.gz \$BASE_DIR/R/project_config.R hamlets_only.pop.txt 6
		"""
}
/* 2) Admixture section ============== */

/* 2a) Admixture (global) -------------- */
process plink12 {
 label 'L_20g2h_plink12'

 input:
 set vcfId, file( vcf ) from vcf_admx

 output:
 set file( "hapmap.ped" ), file( "hapmap.map" ), file( "hapmap.nosex" ) into admx_plink

 script:
 """
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

process admixture {
    label 'L_78g10h_admixture'
    publishDir "2_analysis/admixture/", mode: 'symlink'

    input:
    set val( x ), file( ped ), file( map ), file( nosex ) from admx_prep

    output:
    set file( "hapmap.${x}.Q" ), file( "hapmap.${x}.P" ) into admx_output
		file( "log${x}.out" ) into admx_log

    script:
    """
    admixture --cv ${ped} ${x} | tee log${x}.out
    """
}

process admixture_log {
  label 'L_loc_admixture_log'
  publishDir "2_analysis/admixture/", mode: 'symlink'

  input:
  file( logs ) from admx_log.collect()

  output:
  file( "admixture_report.txt" ) into admxR_output
  script:
  """
  grep -h CV log*.out > admixture_report.txt
  """
}

/* 2b) Admixture (local) -------------- */
process plink12_loc {
 label 'L_20g2h_plink12_loc'

 input:
 set val( loc ), file( vcf ), file( pop ) from vcf_loc_admix

 output:
 set val( loc ), file( "hapmap.${loc}.ped" ), file( "hapmap.${loc}.map" ), file( "hapmap.${loc}.nosex" ) into admx_loc_plink

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
    publishDir "2_analysis/admixture/${loc}", mode: 'symlink'

    input:
    set val( x ), val( loc ), file( ped ), file( map ), file( nosex ) from admx_loc_prep

    output:
    set val( loc ), file( "*.Q" ), file( "*.P" ) into admx_loc_output
		set val( loc ),file( "*.out" ) into admx_loc_log

    script:
    """
    admixture --cv ${ped} ${x} | tee log${x}.out
    """
}

process admixture_loc_log {
  label 'L_loc_admixture_log_loc'
  publishDir "2_analysis/admixture/${loc}", mode: 'symlink'

  input:
  set val( loc ), file( logs ) from admx_loc_log.collect()

  output:
  file( "admixture_report.${loc}.txt" ) into admxR_output

  script:
  """
  grep -h CV log*.out > admixture_report.${loc}.txt
  """
}
/* ============== */
/*
process admixture_plot {
  label 'L_20g2h_admixture_plot'
  publishDir "figures/admixture/", mode: 'symlink'

  input:
  file( report ) from admxR_output
  file( qs ) from admxQ_output.collect()

  output:
  file( "${name}.admixture.pdf") into admxR_plot

  script:
  """
  Rscript --vanilla \$BASE_DIR/R/plot_admixture.R ${report} \$BASE_DIR/vcf_samples.txt ${name}
  """
}
*/
/* 3) fasttree section ============== */
process vcf2geno {
  label 'L_20g15h_vcf2geno'

  input:
  set vcfId, file( vcf ) from vcf_geno

  output:
  file( "output.geno.gz" ) into geno_output

  script:
  """
  python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
    -i ${vcf[0]} | gzip > output.geno.gz
  """
}

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
}

process fasttree {
  label 'L_105g30h_fasttree'
  publishDir "2_analysis/fasttree/", mode: 'symlink'

  input:
  file( geno ) from snp_gene_tree

  output:
  file( "output.SNP.tree" ) into ( fasttree_output )

  script:
  """
  python \$SFTWR/genomics_general/genoToSeq.py -g ${geno} \
      -s output.SNP.phylip \
      -f phylip \
      --splitPhased
  fasttree -nt output.SNP.phylip > output.SNP.tree
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