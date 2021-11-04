#!/usr/bin/env nextflow

// git 13.1
// Open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_snps_ch; vcf_snps_ch2 }

// git 13.2
// Open the allBP data set (will be expanded x 24 LGs)
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.non_ref.vcf.{gz,gz.tbi}")
	.set{ vcf_allbp_ch }

// git 13.3
Channel
	.from( [ 1, 5, 10, 50 ] )
	.set{ window_size_ch }

// git 13.4
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

// git 13.5
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

// git 13.6
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

// git 13.7
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ lg_ch }

// git 13.8
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

// git 13.9
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

// git 13.10
// Compile summary table
process complie_window_stats {
	label 'L_20g2h_windows_stats'
	publishDir "../../2_analysis/window_stats/window_stats/", mode: 'copy' 

	input:
	set val( kb_size ), file( windows ) from coverage_ch.concat( coverage_allbp_ch ).groupTuple()

	output:
	set val( kb_size ), file( "window_stats.${kb_size}kb.tsv.gz" ) into window_out_ch

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

// ----------------------- DISCLAIMER ----------------------
// form here on, this pipeline was not actually run using
//  nextflow, but managed manually
// ---------------------------------------------------------

// git 13.11
// Subset to 5000 random windows
process subset_windows {
	label 'L_loc_subset_windows'

	input:
	set val( kb_size ), file( windows ) from window_out_ch.filter({ it[1] == 5 })

	output:
	set val( kb_size ), file( "5000x_${kb_size}kb_v1.bed" ) into window_subset_ch

	script:
	"""
	#!/usr/bin/env Rscript
	library(hypogen)   # attaches tidyverse automatically

	### Draw windows with coverage (sequence contiguity) and SNP-based cutoffs

	stats5 <- read_tsv("window_stats.${kb_size}kb.tsv.gz")

	set.seed(64)

	selected_windows <- stats5 %>%
	filter(REL_COV >= 0.80 &
			COV_HYP >= 50) %>%
	sample_n(5000) %>%
	arrange(CHROM, START)

	write_tsv(selected_windows %>% select(CHROM, START, END), file = "5000x_${kb_size}kb_v1.bed")
	"""
}

// git 13.12
// index for sub-channels
Channel.from( 1..50 ).into{ loop50_idx_ch1; loop50_idx_ch2; loop50_idx_ch3; loop50_idx_ch4 }

// git 13.13
// Extract windows from all BP
process extract_windows {
	label 'L_20g2h_extract_windows'

	input:
	set val( idx ), file( windows ) from loop50_idx_ch1.combine( window_out_ch )

	output:
	file( "*_v1_all.vcf.gz" ) into window_extract_loop
	
	script:
	"""
	PER_TASK=100

	START_NUM=\$(( (${idx} - 1) * \$PER_TASK + 1 ))
	END_NUM=\$(( ${idx} * \$PER_TASK ))

	echo This is task ${idx}, which will do runs \$START_NUM to \$END_NUM

	for (( run=\$START_NUM ; run<=END_NUM ; run++ ))
		do
		echo This is task ${idx}, run number \$run
		chr=\$(awk -v line="\$run" 'BEGIN { FS = "\\t" } ; NR==line+1 { print \$1 }' 5000x_5kb_v1.bed)
		sta=\$(awk -v line="\$run" 'BEGIN { FS = "\\t" } ; NR==line+1 { print \$2 }' 5000x_5kb_v1.bed)
		end=\$(awk -v line="\$run" 'BEGIN { FS = "\\t" } ; NR==line+1 { print \$3 }' 5000x_5kb_v1.bed)
		printf -v i "%04d" \$run

		vcftools \
			--gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/byLG/filterd.allBP."\$chr".vcf.gz \
			--chr "\$chr" \
			--from-bp "\$sta" \
			--to-bp "\$end" \
			--recode --stdout | \
			grep -v "##" | bgzip > window_"\$i"_v1_all.vcf.gz
	done
	"""
}

// git 13.14
// index for sub-channels
Channel.from( 1..20 ).set{ loop20_idx_ch }
// bundle the distributed sub-channels
window_extract_loop.collect().map{ [ it ] }.set{ window_extract_ch1, window_extract_ch2 }

// git 13.15
// Remove samples (all / noS data sets)
// !!!! should "_noS.vcf.gz" be "v1_noS.vcf.gz" ?
process remove_samples {
	label 'L_20g2h_remove_samples'

	input:
	set val( idx ), file( vcf ) from loop20_idx_ch.combine( window_extract_ch1 )
	
	output:
	file( "*_noS.vcf.gz" ) into samples_removed_loop
	
	script:
	"""
	PER_TASK=250

	START_NUM=\$(( (${idx} - 1) * \$PER_TASK + 1 ))
	END_NUM=\$(( ${idx} * \$PER_TASK ))

	echo This is task ${idx}, which will do runs \$START_NUM to \$END_NUM

	for (( run=\$START_NUM ; run<=END_NUM ; run++ ))
		do

		echo This is task ${idx}, run number \$run
		
		printf -v i "%04d" \$run

		vcftools \
			--gzvcf window_"\$i"*.vcf.gz \
			--remove samples_Serr.txt \
			--recode \
			--stdout | bgzip > window_"\$i"_noS.vcf.gz
	done
	"""
}
// bundle the distributed sub-channels
samples_removed_loop.collect().map{ [ it ] }.set{ samples_removed_ch }

// git 13.16
// Convert to Fasta (IUPAC-encoded)
process fasta_convert {
	label 'L_20g2h_convert_to_fasta'

	input:
	set val( idx ), file( vcf_noS ), file( vcf_all ) from loop50_idx_ch2.combine( samples_removed_ch ).combine( window_extract_ch2 )
	
	output:
	file( "*.fas" ) into fasta_convert_loop
	
	script:
	"""
	PER_TASK=100

	START_NUM=\$(( (${idx} - 1) * \$PER_TASK + 1 ))
	END_NUM=\$(( ${idx} * \$PER_TASK ))

	echo This is task ${idx}, which will do runs \$START_NUM to \$END_NUM

	for (( run=\$START_NUM ; run<=END_NUM ; run++ ))
		do
		echo This is task ${idx}, run number \$run
		printf -v i "%04d" \$run
		for j in all noS
			do
				bgzip -cd window_"\$i"_v1_"\$j".vcf.gz | vcf-to-tab > window_"\$i"_v1_"\$j".tab

				perl \$SFTWR/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i window_"\$i"_v1_"\$j".tab > window_"\$i"_v1_"\$j".fas

				rm window_"\$i"_v1_"\$j".tab*
		done
	done
	"""
}
// bundle the distributed sub-channels
fasta_convert_loop.collect().map{ [ it ] }.set{ fasta_convert_ch }

// git 13.17
// Align sequences in windows
process fasta_align {
	label 'L_20g2h_fasta_align'

	input:
	set val( idx ), file( fa ) from loop50_idx_ch3.combine( fasta_convert_ch )
	
	output:
	file( "*.aln" ) into fasta_align_loop
	
	script:
	"""
	PER_TASK=100

	START_NUM=\$(( (${idx} - 1) * \$PER_TASK + 1 ))
	END_NUM=\$(( ${idx} * \$PER_TASK ))

	echo This is task ${idx}, which will do runs \$START_NUM to \$END_NUM

	for (( run=${idx} ; run<=END_NUM ; run++ ))
		do
		echo This is task ${idx}, run number \$run
		printf -v i "%04d" \$run
		for j in all noS
			do
				mafft --auto window_"\$i"_v1_"\$j".fas > "window_"\$i"_v1_"\$j".aln
		done
	done
	"""
}
// bundle the distributed sub-channels
fasta_align_loop.collect().map{ [ it ] }.set{ fasta_align_ch }

// git 13.18
// Infer local trees
process local_trees {
	label 'L_20g2h_local_trees'

	input:
	set val( idx ), file( aln ) from loop50_idx_ch4.combine( fasta_align_ch )
	
	output:
	file( "*.treefile" ) into local_trees_loop
	
	script:
	"""
	PER_TASK=100

	START_NUM=$(( (${idx}- 1) * \$PER_TASK + 1 ))
	END_NUM=$(( ${idx} * \$PER_TASK ))

	echo This is task ${idx}, which will do runs \$START_NUM to \$END_NUM

	for (( run=${idx} ; run<=END_NUM ; run++ ))
		do
		echo This is task ${idx}, run number \$run
		printf -v i "%04d" \$run
		for j in all noS
			do
				iqtree2 -s window_"\$i"_v1_"\$j".aln  --prefix locus_"\$i"_v1_"\$j" -o PL17_160floflo -T 1
		done
	done
	"""
}
// bundle the distributed sub-channels
local_trees_loop.collect().map{ [ it ] }.set{ local_trees_ch }

// git 13.19
// Calculate summary tree
process summary_trees {
	label 'L_20g2h_summary_trees'
	publishDir "../../2_analysis/astral/", mode: 'copy' 

	input:
	file( tree ) from local_trees_ch
	
	output:
	set file( "astral*.tre" ), file( "astral*.log" ) into summary_trees_ch
	
	script:
	"""
	cat ./*_noS.treefile > genetrees_5000x_5kb_v1_noS.tre
	java -jar \$SFTWR/ASTRAL_5.7.5/astral.5.7.5.jar -i genetrees_5000x_5kb_v1_noS.tre -o astral_5000x_5kb_v1_noS.tre 2> astral_5000x_5kb_v1_noS.log

	cat ./*_all.treefile > genetrees_5000x_5kb_v1_all.tre
	java -jar \$SFTWR/ASTRAL_5.7.5/astral.5.7.5.jar -i genetrees_5000x_5kb_v1_all.tre -o astral_5000x_5kb_v1_all.tre 2> astral_5000x_5kb_v1_all.log
	"""
}