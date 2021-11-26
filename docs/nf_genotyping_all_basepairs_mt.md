---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 18) Genotyping III (all callable sites for mtDNA and unplaced contigs)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/18_genotyping_all_basepairs_mt
source ../../sh/nextflow_alias.sh
nf_run_allbp_mt1
```

## Summary

The genotyping procedure is controlled by the [**nextflow**](https://www.nextflow.io/) script `genotyping_all_basepairs_mt.nf` (located under `$BASE_DIR/nf/18_genotyping_all_basepairs_mt`).
Based on an intermediate step from `genotyping.nf` ([git 1.10](git-1-genotyping-i-snps-only.html)), this script produces a data set that includes _all callable sites_  - that is SNPs as well a invariant sites that are covered by sequence (for mtDNA and unplaced contigs).

The genotypes produced by this script are then used in the Serraninae phylogeny.

## Details of `genotyping_all_basepairs_mt.nf`

### Data preparation

The nextflow script starts with a small header and then imports the joint genotyping likelihoods for all samples produced by `genotyping.nf`.

Furthermore a channel is created to call mtDNA and unplaced contigs seperately.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 18.1</span>
<span class="hl slc">// open genotype likelyhoods</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/1_gvcfs/cohort.g.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_cohort <span class="hl opt">}</span>

<span class="hl kwa">Channel</span>
	.from<span class="hl opt">([</span><span class="hl str">&quot;LG_M&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;unplaced&quot;</span><span class="hl opt">])</span>
	.set<span class="hl opt">{</span> lg_mode <span class="hl opt">}</span>
</code>
</pre>
</div>

The samples are jointly genotyped, independently for mtDNA and unplaced contigs and including invariant sites.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 18.2</span>
<span class="hl slc">// actual genotyping step (including invariant sites)</span>
<span class="hl kwa">process</span> joint_genotype_snps <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_88g48h_LGs_genotype&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../1_genotyping/2_raw_vcfs/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> mode <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_cohort.combine<span class="hl opt">(</span> lg_mode <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_sites.</span><span class="hl ipl">${mode}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_sites.</span><span class="hl ipl">${mode}</span><span class="hl str">.vcf.gz.tbi&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> mode <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> all_bp_non_lg_1<span class="hl opt">,</span> all_bp_non_lg_2 <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	if<span class="hl opt">(</span> mode <span class="hl opt">==</span> <span class="hl str">&#39;unplaced&#39;</span> <span class="hl opt">)</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	gatk --java-options &quot;-Xmx85g&quot; \</span>
<span class="hl str">		GenotypeGVCFs \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-XL=LG01 \</span>
<span class="hl str">		-XL=LG02 \</span>
<span class="hl str">		-XL=LG03 \</span>
<span class="hl str">		-XL=LG04 \</span>
<span class="hl str">		-XL=LG05 \</span>
<span class="hl str">		-XL=LG06 \</span>
<span class="hl str">		-XL=LG07 \</span>
<span class="hl str">		-XL=LG08 \</span>
<span class="hl str">		-XL=LG09 \</span>
<span class="hl str">		-XL=LG10 \</span>
<span class="hl str">		-XL=LG11 \</span>
<span class="hl str">		-XL=LG12 \</span>
<span class="hl str">		-XL=LG13 \</span>
<span class="hl str">		-XL=LG14 \</span>
<span class="hl str">		-XL=LG15 \</span>
<span class="hl str">		-XL=LG16 \</span>
<span class="hl str">		-XL=LG17 \</span>
<span class="hl str">		-XL=LG18 \</span>
<span class="hl str">		-XL=LG19 \</span>
<span class="hl str">		-XL=LG20 \</span>
<span class="hl str">		-XL=LG21 \</span>
<span class="hl str">		-XL=LG22 \</span>
<span class="hl str">		-XL=LG23 \</span>
<span class="hl str">		-XL=LG24 \</span>
<span class="hl str">		-XL=LG_M \</span>
<span class="hl str">		-V=</span><span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		-O=intermediate.vcf.gz \</span>
<span class="hl str">		--include-non-variant-sites=true \</span>
<span class="hl str">		--allow-old-rms-mapping-quality-annotation-data</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx85G&quot; \</span>
<span class="hl str">		SelectVariants \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-V=intermediate.vcf.gz \</span>
<span class="hl str">		--select-type-to-exclude=INDEL \</span>
<span class="hl str">		-O=all_sites.</span><span class="hl ipl">${mode}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	rm intermediate.*</span>
<span class="hl str">	&quot;&quot;&quot;</span>
	else if<span class="hl opt">(</span> mode <span class="hl opt">==</span> <span class="hl str">&#39;LG_M&#39;</span> <span class="hl opt">)</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	gatk --java-options &quot;-Xmx85g&quot; \</span>
<span class="hl str">		GenotypeGVCFs \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-L=</span><span class="hl ipl">${mode}</span> <span class="hl str">\</span>
<span class="hl str">		-V=</span><span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		-O=intermediate.vcf.gz \</span>
<span class="hl str">		--include-non-variant-sites=true \</span>
<span class="hl str">		--allow-old-rms-mapping-quality-annotation-data</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx85G&quot; \</span>
<span class="hl str">		SelectVariants \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-V=intermediate.vcf.gz \</span>
<span class="hl str">		--select-type-to-exclude=INDEL \</span>
<span class="hl str">		-O=all_sites.</span><span class="hl ipl">${mode}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	rm intermediate.*</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
</div>

---
