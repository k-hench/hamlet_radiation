#!/usr/bin/env nextflow

// Open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_snps_ch }

// Open the allBP data set (will be expanded x 24 LGs)
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG*.vcf.{gz,gz.tbi}")
	.set{ vcf_allbp_ch }

// Open the pre-defined window positions
Channel
	.fromPath("../../ressources/windows_1kb.bed.gz")
	.set{ windows_ch }

// Subset ALL vcf files (also allBP) by missingnes (max. 10%)
process filter_vcf_missingnes {
	label "L_20g6h_filter_vcf"
	tag "${vcfId}"

	input:
	set  vcfId, file( vcf ) from vcf_snps_ch.concat( vcf_allbp_ch ).map{ [it[0].minus(".vcf"), it[1]]}
	
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

// Coverage of SNPs vcf for SNPdensity, allBP for Ns
process compute_coverage {
	label "L_20g6h_coverage"
	tag "${vcfId}"
	publishDir "../../2_analysis/revPoMo/coverage", mode: 'copy' 

	input:
	set vcfId, file( vcf ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )
	
	output:
	file( "${vcfId}_cov.tsv.gz" ) into coverage_ch

	script:
	"""
	bedtools coverage \
		-a ${window} \
		-b ${vcf[0]} \
		-counts  > ${vcfId}_cov.tsv
	
	gzip ${vcfId}_cov.tsv
	"""
}

// Compile summary table
process complie_window_stats {
	label "L_20g2h_window_stats"
	publishDir "../../2_analysis/revPoMo/", mode: 'copy' 

	input:
	file( windows ) from coverage_ch.collect()

	output:
	file( "window_stats.tsv.gz" ) into final_ch

	script:
	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	data_SNPs <- read_tsv("phased_mac2_cov.tsv.gz",
							col_names = c("CHROM", "START", "END", "COV_SNP"))

	data_allBPs <- 1:24 %>% 
					str_pad(width = 2, pad = 0) %>%
					str_c("filterd.allBP.LG", ., "".tsv.gz") %>%
					map_dfr(.f = function(file){
							read_tsv(file, 
								col_names = c("CHROM", "START", "END", "COV_ALL") %>%
							filter(COV_ALL > 0 )
							}
						)

	data <- data_SNPs %>%
		left_join(data_allBPs, by = c(CHROM = "CHROM", START = "START", END = "END")) %>%
		mutate(SNP_density = round(COV_SNP/ COV_ALL, 2), 
		REL_COV =  round(COV_ALL/ (END-START), 2))
	
	write_tsv(x = data, file = "window_stats.tsv.gz")
	"""
}