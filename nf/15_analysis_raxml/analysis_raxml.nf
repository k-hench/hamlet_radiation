#!/usr/bin/env nextflow
//
// ----------------------- DISCLAIMER ----------------------
// this pipeline was not actually run using nexflow,
// but manged manually
// ---------------------------------------------------------

// Hamlet Phylogeny
// ---------------

// git 15.1
// Open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_hypo_whg_ch; vcf_serr_whg_ch }

// git 15.2
process hypo_whg_genotypes {
	input:
	set  vcfId, file( vcf ) from vcf_hypo_whg_ch

	output:
	file( "hyp155_n_0.33_mac4_5kb.fas" ) into raxml_hypo_genotypes_ch

	script:
	"""
	# 15.1.1 Remove hybrids and Serranus samples from genotype data (SNPs only)
	vcftools \
	  --gzvcf ${vcf[0]} \
	  --remove samples_155.txt \
	  --recode \
	  --stdout | \
	  gzip > hyp155.vcf.gz

	# 15.1.2 Mask heterozygous genotypes as unknown
	zcat < hyp155.vcf.gz | \
	  sed -e s/"1|0"/".|."/g -e s/"0|1"/".|."/g | \
	  gzip > hyp155_n.vcf.gz

	# 15.1.3 Apply missingness, allele count and distance filters
	vcftools \
	  --gzvcf hyp155_n.vcf.gz \
	  --max-missing 0.33 \
	  --mac 4 \
	  --thin 5000 \
	  --recode \
	  --out hyp155_n_0.33_mac4_5kb

	# 15.1.4 Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl

	vcf-to-tab < hyp155_n_0.33_mac4_5kb.vcf > hyp155_n_0.33_mac4_5kb.tab

	perl ./vcf_tab_to_fasta_alignment.pl -i hyp155_n_0.33_mac4_5kb.tab > hyp155_n_0.33_mac4_5kb.fas
	"""
}

// git 15.3
process hypo_whg_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	file( fas ) from raxml_hypo_genotypes_ch

	output:
	file( "hyp155_n_0.33_mac4_5kb.raxml.support.bs-tbe" ) into raxml_hypo_whg_ch

	script:
	"""
	# 5.1.5 Infer phylogeny
	# Note: number of invariant sites for Felsenstein correction was calculated as number of
	# variant sites in alignment (105,043) / genome-wide proportion of variant sites 
	# (0.05) * genome-wide proportion of invariant sites (0.95)

	raxml-NG --all \
	  --msa ${fas} \
	  --model GTR+G+ASC_FELS{1995817} \
	  --tree pars{20},rand{20} \
	  --bs-trees 100 \
	  --threads 24 \
	  --worker 8 \
	  --seed 123 \
	  --prefix hyp155_n_0.33_mac4_5kb
	"""
}

// RAxML analysis, Serranus-rooted
// ------------------------------
Channel
	.fromPath("../../ressources/samples_hybrids.txt")
	.set{ hybrids_file }

// git 15.5
process serr_whg_genotypes {
	input:
	set  vcfId, file( vcf ), file( hybrids ) from vcf_serr_whg_ch.combine( hybrids_file )

	output:
	file( "hyS_n_0.33_mac4_5kb.fas" ) into raxml_serr_genotypes_ch

	script:
	"""
	# 5.1.1 Remove hybrids from genotype data (SNPs only)
	vcftools \
	  --gzvcf  ${vcf[0]} \
	  --remove ${hybrids} \
	  --recode \
	  --stdout | \
	  gzip > hyS.vcf.gz

	# 5.1.2 Mask heterozygous genotypes as unknown
	zcat < hyS.vcf.gz | \
	  sed -e s/"1|0"/".|."/g -e s/"0|1"/".|."/g | \
	  gzip > hyS_n.vcf.gz

	# 5.1.3 Apply missingness, allele count and distance filters
	vcftools \
	  --gzvcf hyS_n.vcf.gz \
	  --max-missing 0.33 \
	  --mac 4 \
	  --thin 5000 \
	  --recode \
	  --out hyS_n_0.33_mac4_5kb

	# 5.1.4 Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl

	vcf-to-tab < hyS_n_0.33_mac4_5kb.vcf > hyS_n_0.33_mac4_5kb.tab
	
	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i hyS_n_0.33_mac4_5kb.tab > hyS_n_0.33_mac4_5kb.fas
	"""
}

// git 15.6
process serr_whg_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	file( fas ) from raxml_serr_genotypes_ch

	output:
	file( "hyS_n_0.33_mac4_5kb.raxml.support.bs-tbe" ) into raxml_serr_whg_ch

	script:
	"""
	# 5.1.5 Reconstruct phylogeny
	# Note: number of invariant sites for Felsenstein correction was calculated as number of
	# variant sites in alignment (109,660) / genome-wide proportion of variant sites
	# (0.05) * genome-wide proportion of invariant sites (0.95)
	raxml-NG --all \
	  --msa hyS_n_0.33_mac4_5kb.fas \
	  --model GTR+G+ASC_FELS{2083540} \
	  --tree pars{20},rand{20} \
	  --bs-trees 100 \
	  --threads 24 \
	  --worker 4 \
	  --seed 123 \
	  --prefix hyS_n_0.33_mac4_5kb
	"""
}

// Region-specific phylogenies
// ---------------------------

// git 15.7
// bundle allBP files and outlier table
Channel
  .fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG04.vcf.{gz,gz.tbi}")
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .concat(Channel.fromFilePairs("../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}"))
  .merge(Channel.from("LG04_1", "LG12_3", "LG12_4"))
  .combine(Channel.fromPath("../../ressources/focal_outlier.tsv"))
  .set{ vcf_lg_ch }

// git 15.8
// toggle sample modes (with/without serranus outgroup)
Channel.fromPath("../../ressources/samples_155.txt")
  .concat(Channel.fromPath("../../ressources/samples_hybrids.txt"))
  .merge(Channel.from("155", "hyS"))
  .set{ sample_mode_ch }

// git 15.9
process extract_regions {
	input:
	set val( vcfIdx ), file( vcf ), val( outlierId ), file( outlier_file ), file( sample_file ), val( sample_mode ) from vcf_lg_ch.combine( sample_mode_ch )

	output:
	set val( outlierId ), val( sample_mode ), file( "*N.fas" ) into outlier_regions_ch

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
	
	# 6.1.2 Replace unknown character states and asterisks (deletions as encoded by GATK) with "N"
	vcf-to-tab < \${OUT_ALT}_${sample_mode}.vcf | sed -e 's/\\.\\/\\./N\\/N/g' -e 's/[ACGTN\\*]\\/\\*/N\\/N/g' > \${OUT_ALT}_${sample_mode}N.tab

	# 6.1.3 Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl
	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i \${OUT_ALT}_${sample_mode}N.tab > \${OUT_ALT}_${sample_mode}N.fas
	"""
}

// git 15.6
process outlier_regions_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	set val( outlierId ), val( sample_mode ), file( fas ) from outlier_regions_ch

	output:
	file( "*.raxml.support" ) into outlier_results_ch

	script:
	"""
	OUT_ALT=\$(echo ${outlierId} | tr '[:upper:]' '[:lower:]' | sed 's/_/./')

	6.1.4 Reconstruct phylogenies
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

// raxml open questions
// --------------------
// file samples_hybrids.txt needed
// file samples_155.txt needed
// suffix from raml unclear (*.raxml.support vs .raxml.support.bs-tbe)