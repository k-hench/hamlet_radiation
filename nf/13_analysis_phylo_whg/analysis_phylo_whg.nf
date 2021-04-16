#!/usr/bin/env nextflow

// ----------------------- DISCLAIMER ----------------------
// this pipeline was not actually run using nexflow,
// but managed manually
// ---------------------------------------------------------

// Hamlet phylogeny
// ----------------

// git 13.1
// open the SNP data set
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_hypo_whg_ch; vcf_serr_whg_ch }

// RAxML analysis, Serranus-rooted
// -------------------------------
// git 13.2
// open the sample-list (excluding hybrid samples)
Channel
	.fromPath("../../ressources/samples_hybrids.txt")
	.set{ hybrids_file }

// git 13.3
// subset data and convert to fasta for raxml
process serr_whg_genotypes {
	input:
	set vcfId, file( vcf ), file( hybrids ) from vcf_serr_whg_ch.combine( hybrids_file )

	output:
	file( "hyS_n_0.33_mac4_5kb.fas" ) into raxml_serr_genotypes_ch

	script:
	"""
	# Remove hybrids from genotype data (SNPs only)
	vcftools \
	  --gzvcf  ${vcf[0]} \
	  --remove ${hybrids} \
	  --recode \
	  --stdout | \
	  gzip > hyS.vcf.gz

	# Mask heterozygous genotypes as unknown
	zcat < hyS.vcf.gz | \
	  sed -e s/"1|0"/".|."/g -e s/"0|1"/".|."/g | \
	  gzip > hyS_n.vcf.gz

	# Apply missingness, allele count and distance filters
	vcftools \
	  --gzvcf hyS_n.vcf.gz \
	  --max-missing 0.33 \
	  --mac 4 \
	  --thin 5000 \
	  --recode \
	  --out hyS_n_0.33_mac4_5kb

	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl

	vcf-to-tab < hyS_n_0.33_mac4_5kb.vcf > hyS_n_0.33_mac4_5kb.tab
	
	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i hyS_n_0.33_mac4_5kb.tab > hyS_n_0.33_mac4_5kb.fas
	"""
}

// git 13.4
// run raxml (Serranus-rooted)
process serr_whg_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	file( fas ) from raxml_serr_genotypes_ch

	output:
	file( "hyS_n_0.33_mac4_5kb.raxml.support" ) into raxml_serr_whg_ch

	script:
	"""
	# Reconstruct phylogeny
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

// RAxML analysis, floridae-rooted
// -------------------------------
// git 13.5
// open the sample-list (excluding hybrid and Serranus samples)
Channel
	.fromPath("../../ressources/samples_155.txt")
	.set{ hamlet_file }

// git 13.6
// subset data and convert to fasta for raxml
process hypo_whg_genotypes {
	input:
	set vcfId, file( vcf ), file( hamlets ) from vcf_hypo_whg_ch.combine(hamlet_file)

	output:
	file( "hyp155_n_0.33_mac4_5kb.fas" ) into raxml_hypo_genotypes_ch

	script:
	"""
	# Remove hybrid and Serranus samples from genotype data (SNPs only)
	vcftools \
	  --gzvcf ${vcf[0]} \
	  --remove ${hamlets} \
	  --recode \
	  --stdout | \
	  gzip > hyp155.vcf.gz

	# Mask heterozygous genotypes as unknown
	zcat < hyp155.vcf.gz | \
	  sed -e s/"1|0"/".|."/g -e s/"0|1"/".|."/g | \
	  gzip > hyp155_n.vcf.gz

	# Apply missingness, allele count and distance filters
	vcftools \
	  --gzvcf hyp155_n.vcf.gz \
	  --max-missing 0.33 \
	  --mac 4 \
	  --thin 5000 \
	  --recode \
	  --out hyp155_n_0.33_mac4_5kb

	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)
	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl

	vcf-to-tab < hyp155_n_0.33_mac4_5kb.vcf > hyp155_n_0.33_mac4_5kb.tab

	perl ./vcf_tab_to_fasta_alignment.pl -i hyp155_n_0.33_mac4_5kb.tab > hyp155_n_0.33_mac4_5kb.fas
	"""
}

// git 13.7
// run raxml (floridae-rooted)
process hypo_whg_raxml {
	publishDir "../../2_analysis/raxml/", mode: 'copy' 
	
	input:
	file( fas ) from raxml_hypo_genotypes_ch

	output:
	file( "hyp155_n_0.33_mac4_5kb.raxml.support" ) into raxml_hypo_whg_ch

	script:
	"""
	# Infer phylogeny
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
