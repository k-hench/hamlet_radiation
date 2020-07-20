#!/usr/bin/env nextflow
// git 11.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

Channel
	.from( "hamlet_only" , "all" )
	.set{ sample_modes }

// git 11.2
// load focal outlier regions
Channel
	.fromPath("../../ressources/focal_outlier.tsv")
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ chrom:row.chrom, start:row.start, end:row.end, gid:row.gid ] }
	.combine( sample_modes )
	.combine( vcf_ch )
	.set{ starter_ch }

// git 11.3
// subset the genotypes by location
process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
//	set val( outlier ), val( sample_mode ) from starter_ch
	set val( outlier ), val( sample_mode ), val( vcfidx ), file( vcf ) from starter_ch

	output:
	set val( outlier.gid ), val( sample_mode ), file( "${outlier.gid}.${sample_mode}.vcf.gz" ) into ( vcf_filtered )

	script:
	"""
	echo -e "CHROM\\tSTART\\tEND" > outl.bed
	echo -e "${outlier.chrom}\\t${outlier.start}\\t${outlier.end}" >> outl.bed

	echo -e "20478tabhon\\n28393torpan\\ns_tort_3torpan" > outgr.pop

	# check if outgroups need to be dropped
	if [ "${sample_mode}" == "all" ];then
		INDS=""
	else
		INDS="--remove outgr.pop"
	fi

	#vcftools --gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/filterd.allBP.${outlier.chrom}.vcf.gz \
	vcftools --gzvcf ${vcf[0]} \
	  \$INDS \
		--bed outl.bed \
		--recode \
		--stdout | gzip > ${outlier.gid}.${sample_mode}.vcf.gz
	"""
}

// git 11.4
// subset genotypes by LG
process vcf2geno_loc {
	label 'L_2g15m_vcf2geno'

	input:
	set val( gid ), val( sample_mode ), file( vcf ) from vcf_filtered

	output:
	set val( gid ), val( sample_mode ), file( "${gid}.${sample_mode}.geno.gz" ) into ( geno_filtered )

	script:
	"""
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
	  -i ${vcf}| gzip > ${gid}.${sample_mode}.geno.gz
	"""
}

// git 11.17
// create the phylogenies along the sliding window
process twisst_prep {
  label "L_30g2h4x_subset_vcf_whg"
  publishDir "../../2_analysis/sliding_phylo/positions/${loc}/", mode: 'copy'

  input:
  set val( gid ), val( sample_mode ), file( geno ) from ( geno_filtered )

	output:
	set file( "*.trees.gz" ), file( "*.data.tsv" ) into twisst_prep_ch

  script:
   """
		python \$SFTWR/genomics_general/phylo/phyml_sliding_windows.py \
      -g ${geno} \
      --windType sites \
      -S 10000 \
      --prefix ${gid}.${sample_mode}.phyml_bionj \
      --model GTR \
      --optimise n \
		--threads 2
	 """
}
