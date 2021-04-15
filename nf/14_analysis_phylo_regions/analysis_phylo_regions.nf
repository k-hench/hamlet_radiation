#!/usr/bin/env nextflow

// Region-specific phylogenies
// ---------------------------

// git 14.1
// bundle allBP files and outlier table
Channel
  .fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG04.vcf.{gz,gz.tbi}")
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .merge(Channel.from("LG04_1", "LG12_3", "LG12_4"))
  .combine(Channel.fromPath("../../ressources/focal_outlier.tsv"))
  .set{ vcf_lg_ch }

// git 14.2
// toggle sample modes (with/ without Serranus outgroup)
Channel.fromPath("../../ressources/samples_155.txt")
  .concat(Channel.fromPath("../../ressources/samples_hybrids.txt"))
  .merge(Channel.from("155", "hyS"))
  .set{ sample_mode_ch }

// git 14.3
process extract_regions {

	input:
	set val( vcfIdx ), file( vcf ), val( outlierId ), file( outlier_file ), file( sample_file ), val( sample_mode ) from vcf_lg_ch.combine( sample_mode_ch )

	output:
	set file( "*_${sample_mode}.vcf" ), val( outlierId ), val( sample_mode ) into ( vcf_raxml_ch, vcf_pomo_ch )

	script:
	"""
	# Extract regions of interest from genotype data (allBP),
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
	"""
}

// git 14.4
process run_pomo {
	publishDir "../../2_analysis/revPoMo/outlier_regions/", mode: 'copy' 
	
	input:
	set file( vcf ), val( outlierId ), val( sample_mode ) from vcf_raxml_ch

	output:
	file( "*_pop.cf.treefile" ) into pomo_results_ch

	script:
	"""
	OUT_ALT=\$(echo ${outlierId} | tr '[:upper:]' '[:lower:]' | sed 's/_/./')

	# Convert to fasta format (Python scripts available at https://github.com/simonhmartin/genomics_general), picked up from 6.1.1 output
	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py -i ${vcf} > \${OUT_ALT}_${sample_mode}.geno

	python \$SFTWR/genomics_general/genoToSeq.py \
		-g \${OUT_ALT}_${sample_mode}.geno \
		-s \${OUT_ALT}_${sample_mode}.fas \
		-f fasta \
		--splitPhased
	
	# Reformat sample ids to provide population prefixes for cflib
	sed -e 's/-/_/g' -e 's/>\(.*\)\([a-z]\{6\}\)_\([AB]\)/>\2-\1_\3/g' \${OUT_ALT}_${sample_mode}.fas > \${OUT_ALT}_${sample_mode}_p.fas

	# Convert to allele frequency format (cflib library available at https://github.com/pomo-dev/cflib)
	\$SFTWR/cflib/FastaToCounts.py \${OUT_ALT}_${sample_mode}_p.fas \${OUT_ALT}_${sample_mode}_pop.cf

	# IQTREE analysis under PoMo model
	iqtree2 \
		-nt 16 \
		-s \${OUT_ALT}_${sample_mode}_pop.cf \
		-m HKY+F+P+N9+G4 \
		-b 100 \
		--tbe
	"""
}

// Region-specific phylogenies
// ---------------------------
// git 14.5
process conversion_raxml {
	input:
	set file( vcf ), val( outlierId ), val( sample_mode ) from vcf_pomo_ch
	
	output:
	set val( outlierId ), val( sample_mode ), file( "*N.fas" ) into outlier_regions_ch

	script:
	"""
	OUT_ALT=\$(echo ${outlierId} | tr '[:upper:]' '[:lower:]' | sed 's/_/./')

	# Replace unknown character states and asterisks (deletions as encoded by GATK) with "N"
	vcf-to-tab < ${vcf} | sed -e 's/\\.\\/\\./N\\/N/g' -e 's/[ACGTN\\*]\\/\\*/N\\/N/g' > \${OUT_ALT}_${sample_mode}N.tab

	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl
	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i \${OUT_ALT}_${sample_mode}N.tab > \${OUT_ALT}_${sample_mode}N.fas
	"""
}

// git 14.6
process run_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	set val( outlierId ), val( sample_mode ), file( fas ) from outlier_regions_ch

	output:
	file( "*.raxml.support" ) into outlier_results_ch

	script:
	"""
	OUT_ALT=\$(echo ${outlierId} | tr '[:upper:]' '[:lower:]' | sed 's/_/./')

	# Reconstruct phylogenies
	raxml-NG --all \
		--msa ${fas} \
		--model GTR+G \
		--tree pars{10},rand{10} \
		--bs-trees 100 \
		--threads 24 \
		--worker 8 \
		--seed 123 \
		--prefix \${OUT_ALT}_${sample_mode}N
	"""
}