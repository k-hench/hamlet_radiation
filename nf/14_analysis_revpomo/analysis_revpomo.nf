#!/usr/bin/env nextflow
// git 14.1
// Open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_snps_ch }

// git 14.2
// Prepare LG channel for vcf subsetting
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24'))
	.map{ "LG" + it }
	.set{ lg_ch }

// git 14.3
// Subset snp data set by LG
process subset_snps_by_lg {
	label "L_20g2h_subset_lg"
	tag "${vcfId}"

	input:
	set  vcfId, file( vcf ), val ( lg ) from vcf_snps_ch.map{ [it[0].minus(".vcf"), it[1]]}.combine( lg_ch )
	
	output:
	set val( "${vcfId}.${lg}.vcf" ), file( "${vcfId}.${lg}.vcf.gz*" ) into vcf_snps_lg_ch

	script:
	"""
	vcftools \
		--gzvcf ${vcf[0]} \
		--chr ${lg} \
		--recode \
		--stdout | \
		bgzip > ${vcfId}.${lg}.vcf.gz
	
	tabix ${vcfId}.${lg}.vcf.gz
	"""
}

// git 14.4
// Open the allBP data set (will be expanded x 24 LGs)
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG*.vcf.{gz,gz.tbi}")
	.set{ vcf_allbp_ch }

// git 14.5
// Open the pre-defined window positions
Channel
	.fromPath("../../ressources/windows_1kb.bed.gz")
	.set{ windows_ch }

// git 14.6
// Subset ALL vcf files (also allBP) by missingnes (max. 10%)
process filter_vcf_missingnes {
	label "L_20g6h_filter_vcf"
	tag "${vcfId}"

	input:
	set  vcfId, file( vcf ) from vcf_snps_lg_ch.concat( vcf_allbp_ch ).map{ [it[0].minus(".vcf"), it[1]]}
	
	output:
	set val( vcfId ),  file( "*_single_ind.vcf.gz*" ) into vcf_snps_filterd_ch

	script:
	"""
	vcftools \
		--gzvcf ${vcf[0]} \
		--max-missing 0.9 \
		--recode \
		--stdout | \
		sed "s/<NON_REF>/./g" | \
		bgzip > ${vcfId}_filtered.vcf.gz
	
	tabix ${vcfId}_filtered.vcf.gz

	vcfsamplenames ${vcfId}_filtered.vcf.gz | head -n 1 > first_ind.txt

	vcftools \
		--gzvcf ${vcfId}_filtered.vcf.gz \
		--keep first_ind.txt \
		--recode \
		--stdout | \
		bgzip > ${vcfId}_single_ind.vcf.gz
	
	tabix ${vcfId}_single_ind.vcf.gz
	"""
}

// git 14.7
// Coverage of SNPs vcf for SNPdensity, allBP for Ns
process compute_coverage {
	label "L_50g3h_coverage"
	tag "${vcfId}"
	publishDir "../../2_analysis/revPoMo/coverage", mode: 'copy' 

	input:
	set vcfId, file( vcf ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )
	
	output:
	file( "${vcfId}_cov.tsv.gz" ) into coverage_ch

	script:
	"""
	LG=\$( echo ${vcfId} | sed 's/.*\\(LG[0-9]\\{2\\}\\)/\\1/' )
	echo -e "CHROM\\tSTART\\tEND" > windows_1kb.\$LG.bed

	zcat ${window} | \
		grep \$LG >> windows_1kb.\$LG.bed
	
	gzip windows_1kb.\$LG.bed

	bedtools coverage \
		-a windows_1kb.\$LG.bed.gz \
		-b ${vcf[0]} \
		-counts  > ${vcfId}_cov.tsv
	
	gzip ${vcfId}_cov.tsv
	"""
}

// git 14.8
// Compile summary table
process compile_window_stats {
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

	data_SNPs <- 1:24 %>% 
					str_pad(width = 2, pad = 0) %>%
					str_c("phased_mac2.LG", ., "_cov.tsv.gz") %>%
					map_dfr(.f = function(file){
							read_tsv(file, 
								col_names = c("CHROM", "START", "END", "COV_SNP")) %>%
							filter(COV_SNP > 0 )
								} ) 
					
	data_allBPs <- 1:24 %>% 
					str_pad(width = 2, pad = 0) %>%
					str_c("filterd.allBP.LG", ., "_cov.tsv.gz") %>%
					map_dfr(.f = function(file){
							read_tsv(file, 
								col_names = c("CHROM", "START", "END", "COV_ALL")) %>%
							filter(COV_ALL > 0 )
							}
						)

	data <- data_SNPs %>%
		left_join(data_allBPs, by = c(CHROM = "CHROM", START = "START", END = "END")) %>%
		mutate(SNP_density = COV_SNP/ COV_ALL, 
				REL_COV =  COV_ALL/ (END-START))
	
	data %>% write_tsv("window_stats.tsv.gz")
	"""
}

// ----------------------- DISCLAIMER ----------------------
// from this point on, the following steps were not actually
// run using nexflow, but manged manually
// ---------------------------------------------------------

// git 14.9
Channel
	.from{ 1, 2 }
	.set{ revpomo_run_ch }

// git 14.10
process select_random_windows {
	input:
	set file( window_stats ), val( revpomo_run ) from final_ch.combine( revpomo_run_ch )

	output:
	set val( revpomo_run ), file( "5k1_80-80_${revpomo_run}.bed" ) into selected_windows_ch

	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/snp_density_ci.R ${window_stats} 5k1_80-80_${revpomo_run}.bed
	"""
}

// git 14.11
Channel
	.fromPath("../../ressources/samples_155.txt")
	.set{ 155file_ch }

// git 14.12
process extract_windows_from_genotypes {
	publishDir "../../2_analysis/revPoMo/", mode: 'copy' 
	
	input:
	set val( revpomo_run ), file( bed ), file( 155file ) from selected_windows_ch.combine( 155file_ch )

	// output:

	script:
	"""
	for i in {01..24} ; do 
		vcftools \
			--gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG\$i.vcf.gz \
			--bed ${bed} \
			--remove-indels \
			--recode \
			--out 5k1_80-80_${revpomo_run}.LG\$i;
	done

	## ---------------- ##
	#  vcfs now merged?  #
	## ---------------- ##

	vcftools \
		--vcf 5k1_80-80_1h.vcf \
		--remove ${155file} \
		--recode \
		--out 5k1_80-80_${revpomo_run}h_155

	# 5.3.4 Convert to fasta format (Python scripts available at https://github.com/simonhmartin/genomics_general)
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py -i 5k1_80-80_${revpomo_run}h_155.vcf > 5k1_80-80_${revpomo_run}h_155.geno

	python \$SFTWR/genomics_general/genoToSeq.py \
		-g 5k1_80-80_${revpomo_run}h_155.geno \
		-s 5k1_80-80_${revpomo_run}h_155.fas \
		-f fasta \
		--splitPhased
	
	# 5.3.5 Reformat sample ids to provide population prefixes for cflib
	sed -e 's/-/_/g' -e 's/>\(.*\)\([a-z]\{6\}\)_\([AB]\)/>\2-\1_\3/g' 5k1_80-80_${revpomo_run}h_155.fas > 5k1_80-80_${revpomo_run}h_155p.fas

	# 5.3.6 Convert to allele frequency format
	# (cflib library available at https://github.com/pomo-dev/cflib)
	\$SFTWR/cflib/FastaToCounts.py 5k1_80-80_${revpomo_run}h_155p.fas 5k1_80-80_${revpomo_run}h_155pop.cf

	# 5.3.7 IQTREE analysis under PoMo model
	iqtree2 \
		-T AUTO \
		-s 5k1_80-80_${revpomo_run}h_155pop.cf \
		-m HKY+F+P+N9+G4 \
		-allnni \
		-b 100

	### I have never seen these output files I think
	"""
}

// Region-specific phylogenies
// ---------------------------

// git 14.13
// bundle allBP files and outlier table
Channel
  .fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG04.vcf.{gz,gz.tbi}")
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .merge(Channel.from("LG04_1", "LG12_3", "LG12_4"))
  .combine(Channel.fromPath("../../ressources/focal_outlier.tsv"))
  .set{ vcf_lg_ch }

// git 14.14
// toggle sample modes (with/without serranus outgroup)
Channel.fromPath("../../ressources/samples_155.txt")
  .concat(Channel.fromPath("../../ressources/samples_hybrids.txt"))
  .merge(Channel.from("155", "hyS"))
  .set{ sample_mode_ch }

// git 14.15
process extract_regions {
	publishDir "../../2_analysis/revPoMo/outlier_regions/", mode: 'copy' 
	
	input:
	set val( vcfIdx ), file( vcf ), val( outlierId ), file( outlier_file ), file( sample_file ), val( sample_mode ) from vcf_lg_ch.combine( sample_mode_ch )

	output:
	file( "*_pop.cf.treefile" ) into pomo_results_ch

	script:
	"""
	# 6.1.1 Extract regions of interest from genotype data (allBP),
	# remove hybrid / Serranus samples and indels; simplify headers

	head -n 1 ${outlier_file} | cut -f 1-3 > outlier.bed
	grep ${outlierId} ${outlier_file} | cut -f 1-3 >> outlier.bed

	OUT_ALT=\$(echo ${outlierId} | tr '[:upper:]' '[:lower:]' | sed 's/_/./')

	vcftools --gzvcf \
	  ${vcf[0]} \
	  --bed outlier.bed \
	  --remove-indels \
	  --remove ${sample_file} \
	  --recode \
	  --stdout | \
	  grep -v '##' > \${OUT_ALT}_${sample_mode}.vcf   # 113,099 / 146,663 / 97,653 sites

	# 6.3.1 Convert to fasta format (Python scripts available at https://github.com/simonhmartin/genomics_general), picked up from 6.1.1 output
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py -i \${OUT_ALT}_${sample_mode}.vcf > \${OUT_ALT}_${sample_mode}.geno

	python \$SFTWR/genomics_general/genoToSeq.py \
		-g \${OUT_ALT}_${sample_mode}.geno \
		-s \${OUT_ALT}_${sample_mode}.fas \
		-f fasta \
		--splitPhased
	
	# 6.3.2 Reformat sample ids to provide population prefixes for cflib
	sed -e 's/-/_/g' -e 's/>\(.*\)\([a-z]\{6\}\)_\([AB]\)/>\2-\1_\3/g' \${OUT_ALT}_${sample_mode}.fas > \${OUT_ALT}_${sample_mode}_p.fas

	# 6.3.3 Convert to allele frequency format (cflib library available at https://github.com/pomo-dev/cflib)
	\$SFTWR/cflib/FastaToCounts.py \${OUT_ALT}_${sample_mode}_p.fas \${OUT_ALT}_${sample_mode}_pop.cf

	# 6.3.4 IQTREE analysis under PoMo model
	iqtree2 \
		-nt 16 \
		-s \${OUT_ALT}_${sample_mode}_pop.cf \
		-m HKY+F+P+N9+G4 \
		-b 100 \
		--tbe
	"""
}

// revpomo open questions
// --------------------
// file snp_density_ed.R needed
// what does "VCF headers were simplified" mean (# 5.3.4)
// revpomo whg never seen