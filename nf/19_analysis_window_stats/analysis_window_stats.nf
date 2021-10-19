#!/usr/bin/env nextflow

// Open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_snps_ch }

// Open the allBP data set (will be expanded x 24 LGs)
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.non_ref.vcf.{gz,gz.tbi}")
	.set{ vcf_allbp_ch }

Channel
	.from( [ 1, 5, 10, 50 ] )
	.set{ window_size_ch }


// Compile summary table
process segment_windows {
	label 'L_loc_slice_windows'
	publishDir "../../2_analysis/window_stats/windows/", mode: 'copy' 

	input:
	val( kb_size ) from window_size_ch

	output:
	set val( kb_size ), file( "windows_${kb_size}kb.bed.gz" ) into ( windows_ch, windows_ch2, windows_ch3 )

	script:
	"""
	#!/usr/bin/env Rscript
	library(hypogen)

	x_bp <- ${kb_size} * 1000
	
	window_non_overlap <- function(CHROM, LENGTH, windowsize = x_bp ){
	  tibble(CHROM = CHROM, 
	         START = seq(from = 1, to = LENGTH, by = windowsize) - 1, # (produces overlap of one which is needed for bedtools)
	         END = lead(START, default = LENGTH) ) }

	hypo_karyotype %>%
	  select(CHROM, LENGTH) %>%
	  pmap_dfr(window_non_overlap) %>%
	  write_tsv(file = "windows_${kb_size}kb.bed.gz")
	"""
}


/*
// Subset ALL vcf files (also allBP) by missingnes (max. 10%)
process filter_vcf_missingnes {
	label 'L_32g12h_filter_missingnes'

	input:
	set  val( vcfId ), val( lg ) from vcfId_allbp_ch.combine( lg_ch )
	
	output:
	set val( "${vcfId}.LG${lg}" ),  file( "*_filtered.vcf.gz*" ) into vcf_snps_filterd_ch

	script:
	"""
	vcftools \
		--gzvcf \BASE_DIR/1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG${lg}.vcf.gz \
		--max-missing 0.9 \
		--recode \
		--stdout | bgzip > ${vcfId}.LG${lg}_filtered.vcf.gz
	
	tabix ${vcfId}.LG${lg}_filtered.vcf.gz
	"""
}*/

// Subset ALL vcf files (also allBP) by missingnes (max. 10%)
/*
process filter_vcf_missingnes {
	label 'L_32g48h_filter_missingnes'

	input:
	set  val( vcfId ), file( vcf ) from vcf_snps_ch.concat( vcf_allbp_ch ).map{ [it[0].minus(".vcf"), it[1]]}
	
	output:
	set val( vcfId ),  file( "*_filtered.vcf.gz*" ) into vcf_snps_filterd_ch

	script:
	"""
	vcftools \
		--gzvcf ${vcf[0]} \
		--max-missing 0.9 \
		--recode \
		--stdout | bgzip > ${vcfId}_filtered.vcf.gz
	
	tabix ${vcfId}_filtered.vcf.gz
	"""
}
*/

process filter_hamlets {
	label 'L_20g2h_outgroup_drop'

	input:
	set  val( vcfId ), file( vcf ) from vcf_snps_ch2

	output:
	set val( "hamlets_only" ), file( "hamlets_only.vcf.gz*" ) into vcf_hamlets_ch

	script:
	"""
	echo -e "20478tabhon\\n28393torpan\\ns_tort_3torpan" > outgroup.pop

	vcftools  --gzvcf ${vcf[0]} \
		--remove outgroup.pop \
		--mac 1 \
		--recode \
		--stdout | \
		bgzip > hamlets_only.vcf.gz
	
	tabix hamlets_only.vcf.gz
	"""
}

// Coverage of SNPs vcf for SNPdensity, allBP for Ns
process compute_coverage {
	label 'L_140g1h_coverage'
	publishDir "../../2_analysis/window_stats/coverages/", mode: 'copy' 

	input:
	//set vcfId, file( vcf ), val( kb_size ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )
	set  val( vcfId ), file( vcf ), val( kb_size ), file( window ) from vcf_snps_ch.concat( vcf_hamlets_ch ).map{ [it[0].minus(".vcf"), it[1]]}.combine( windows_ch )
	
	output:
	set val( kb_size ), file( "${vcfId}.${kb_size}kb_cov.tsv.gz" ) into coverage_ch

	script:
	"""
	bedtools coverage \
		-a ${window} \
		-b ${vcf[0]} \
		-counts | gzip > ${vcfId}.${kb_size}kb_cov.tsv.gz
	"""
}

Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ lg_ch }

process subset_allBP {
	label 'L_140g10h_coverage'
	publishDir "../../1_genotyping/3_gatk_filtered/non_ref_byLG/", mode: 'copy' 

	input:
	//set vcfId, file( vcf ), val( kb_size ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )
	set  val( lg ) from lg_ch

	output:
	set val( lg ), file( "filterd.allBP.non_ref.LG${lg}.vcf.gz" ), file( "filterd.allBP.non_ref.LG${lg}.vcf.gz.tbi" ) into allbp_non_ref_ch

	script:
	"""
	vcftools \
		--gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered//filterd.allBP.non_ref.vcf.gz \
		--chr LG${lg} \
		--recode \
		--stdout | bgzip > filterd.allBP.non_ref.LG${lg}.vcf.gz

	tabix filterd.allBP.non_ref.LG${lg}.vcf.gz
	"""
}


//vcf_snps_ch.concat( vcf_allbp_ch ).map{ [it[0].minus(".vcf"), it[1]]}
process compute_coverage_allBP {
	label 'L_140g1h_coverage'
	publishDir "../../2_analysis/window_stats/coverages/", mode: 'copy' 

	input:
	//set vcfId, file( vcf ), val( kb_size ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )
	set  val( lg ), file( vcf ), file( tbi ), val( kb_size ), file( window ) from allbp_non_ref_ch.combine( windows_ch2 )
	
	output:
	set val( kb_size ), file( "filterd.allBP.LG${lg}.${kb_size}kb_cov.tsv.gz" ) into coverage_allbp_ch

	script:
	"""
	echo -e "CHROM\\tSTART\\tEND" > bed_LG${lg}.bed

	zgrep "LG${lg}" ${window} >> bed_LG${lg}.bed
	gzip bed_LG${lg}.bed

	bedtools coverage \
		-a bed_LG${lg}.bed.gz \
		-b ${vcf} \
		-counts | gzip > filterd.allBP.LG${lg}.${kb_size}kb_cov.tsv.gz
	"""
}


// Compile summary table
process complie_window_stats {
	label 'L_20g2h_windows_stats'
	publishDir "../../2_analysis/window_stats/window_stats/", mode: 'copy' 

	input:
	set val( kb_size ), file( windows ) from coverage_ch.concat( coverage_allbp_ch ).groupTuple()

	output:
	file( "window_stats.${kb_size}kb.tsv.gz" ) into final_ch

	script:
	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	data_SNPs <- read_tsv("phased_mac2.${kb_size}kb_cov.tsv.gz",
						  col_names = c("CHROM", "START", "END", "COV_SNP"))

	data_HYP <- read_tsv("hamlets_only.${kb_size}kb_cov.tsv.gz",
						  col_names = c("CHROM", "START", "END", "COV_HYP"))

	all_bp_files <- dir(pattern = "filterd.allBP.LG*")

	data_allBPs <- map_dfr(all_bp_files, .f = function(x){read_tsv(x, col_names = c("CHROM", "START", "END", "COV_ALL"))})

	data <- data_SNPs %>%
		left_join(data_HYP, by = c(CHROM = "CHROM", START = "START", END = "END"))  %>%
		left_join(data_allBPs, by = c(CHROM = "CHROM", START = "START", END = "END"))  %>%
		filter(COV_ALL > 0 ) %>%
		mutate(SNP_density = round(COV_SNP/ COV_ALL, 2), 
		REL_COV =  round(COV_ALL/ (END-START), 2))
	
	write_tsv(x = data, file = "window_stats.${kb_size}kb.tsv.gz")
	"""
}
