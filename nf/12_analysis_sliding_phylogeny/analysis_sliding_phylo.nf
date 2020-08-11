#!/usr/bin/env nextflow
// git 12.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

Channel
	.from( "hamlet_only" , "all" )
	.set{ sample_modes }

Channel
	.from( 50, 100, 150, 200 )
	.set{ window_ch }

// git 12.2
// load focal outlier regions
Channel
	.fromPath("../../ressources/focal_outlier.tsv")
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ chrom:row.chrom, start:row.start, end:row.end, gid:row.gid ] }
	.combine( sample_modes )
	.combine( vcf_ch )
	.set{ starter_ch }

// git 12.3
// subset the genotypes by location
process subset_vcf_by_location {
	label "L_20g2h_subset_vcf"

	input:
//	set val( outlier ), val( sample_mode ) from starter_ch
	set val( outlier ), val( sample_mode ), val( vcfidx ), file( vcf ) from starter_ch

	output:
	set val( outlier.gid ), val( sample_mode ), file( "${outlier.gid}.${sample_mode}.vcf.gz" ), file( "outl.bed")  into ( vcf_filtered )

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

	# vcftools --gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/filterd.allBP.${outlier.chrom}.vcf.gz

	vcftools --gzvcf ${vcf[0]} \
	  \$INDS \
		--mac 3 \
		--bed outl.bed \
		--recode \
		--stdout | gzip > ${outlier.gid}.${sample_mode}.vcf.gz
	"""
}

// git 12.4
// subset genotypes by LG
process vcf2geno_loc {
	label 'L_2g15m_vcf2geno'

	input:
	set val( gid ), val( sample_mode ), file( vcf ), file( bed )  from vcf_filtered

	output:
	set val( gid ), val( sample_mode ), file( "${gid}.${sample_mode}.geno.gz" ), file( bed ) into ( geno_filtered  )

	script:
	"""
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
	  -i ${vcf}| gzip > ${gid}.${sample_mode}.geno.gz
	"""
}

// git 12.17
geno_filtered
	.combine( window_ch )
	.into{ phyml_input_ch; raxml_input_ch }

// git 12.17
// create the phylogenies along the sliding window
process phyml_slide {
  label "L_30g2h4x_subset_vcf_whg"
  publishDir "../../2_analysis/sliding_phylo/", mode: 'copy'

  input:
  set val( gid ), val( sample_mode ), file( geno ), file( bed ), val( win ) from phyml_input_ch

	output:
	set file( "*.trees.gz" ), file( "*.data.tsv" ) into phyml_out

  script:
   """
	 tail -n 1 ${bed} | \
	   awk '{l = \$1; s = \$2; e = \$3; while (s < e) { print l" "s" "s+999 ; s = s + 1000} }' \
		 > steps.bed

	python \$SFTWR/genomics_general/phylo/phyml_sliding_windows.py \
      -g ${geno} \
			--windType sites \
      -w ${win} \
      --prefix ${gid}.${sample_mode}.${win}w.phyml_bionj \
      --model GTR \
      --optimise n \
		--threads 2
	 """
}

process raxml_slide {
  label "L_30g2h4x_subset_vcf_whg"
  publishDir "../../2_analysis/sliding_phylo/", mode: 'copy'

  input:
  set val( gid ), val( sample_mode ), file( geno ), file( bed ), val( win ) from raxml_input_ch

	output:
	set file( "*.trees.gz" ), file( "*.data.tsv" ) into raxml_out

  script:
   """
	 tail -n 1 ${bed} | \
	   awk '{l = \$1; s = \$2; e = \$3; while (s < e) { print l" "s" "s+999 ; s = s + 1000} }' \
		 > steps.bed

	python \$SFTWR/genomics_general/phylo/raxml_sliding_windows.py \
      -g ${geno} \
			--windType sites \
      -w ${win} \
      --prefix ${gid}.${sample_mode}.${win}w.raxml_bionj \
      --model GTRCAT \
		--threads 2
	 """
}
