---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 17) Analysis XV (dstats)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/17_analysis_dstats
nextflow run analysis_dstats.nf -c ../../nextflow.config -resume
```

## Summary

The <span style="color:red;">...</span> are computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_dstats.nf` (located under `$BASE_DIR/nf/17_analysis_dstats/`).
It takes the <span style="color:red;">...</span> and computes <span style="color:red;">...</span>.
Below is an overview of the steps involved in the analysis.

## Details of `analysis_dstats.nf`

### Setup

The nextflow script starts by <span style="color:red;">...</span>

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow

<span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// this pipeline was not actually run using nextflow,</span>
<span class="hl slc">// but managed manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>

<span class="hl slc">// git 17.1</span>
<span class="hl slc">// load genotypes</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> genotypes_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.2</span>
<span class="hl slc">// drop serranus</span>
<span class="hl kwa">process</span> drop_serranus <span class="hl opt">{</span>
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> genotypes_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp_mac2.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> no_serranus_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep &quot;tor</span><span class="hl esc">\\</span><span class="hl str">|tab</span><span class="hl esc">\\</span><span class="hl str">&quot; &gt; serr.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--remove serr.pop \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt; hyp_mac2.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.4</span>
<span class="hl slc">// LD filtering</span>
<span class="hl kwa">process</span> ld_filter <span class="hl opt">{</span>
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> no_serranus_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp_mac2_ld05.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> ld_filtered_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	bcftools +prune \</span>
<span class="hl str">		-l 0.5 \</span>
<span class="hl str">		-w 50kb \</span>
<span class="hl str"></span>		<span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		-Oz \</span>
<span class="hl str">		-o hyp_mac2_ld05.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.5</span>
<span class="hl slc">// run D trios</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/hyp_sets.txt&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> hyp_sets_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.6</span>
<span class="hl slc">// run D trios</span>
<span class="hl kwa">process</span> run_dtrios <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/dstats/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> sets <span class="hl opt">)</span> <span class="hl kwa">from</span> ld_filtered_ch.combine<span class="hl opt">(</span> hyp_sets_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp_ld05_dtrios_BBAA.txt&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp_ld05_dtrios_Dmin.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> dtrios_results_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	zcat</span> <span class="hl ipl">${vcf}</span> <span class="hl str">&gt; hyp_mac2_ld05.vcf</span>
<span class="hl str"></span>
<span class="hl str">	Dsuite Dtrios \</span>
<span class="hl str">		-c \</span>
<span class="hl str">		-o hyp_ld05_dtrios \</span>
<span class="hl str">		hyp_mac2_ld05.vcf \</span>
<span class="hl str"></span>		<span class="hl ipl">${sets}</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.7</span>
<span class="hl slc">// load predefined species order</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/species_order_alpha.txt&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> spec_order_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 17.8</span>
<span class="hl slc">// Extract significant Dmin, BBAA trios</span>
<span class="hl kwa">process</span> run_correction <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/dstats/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> sets <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> spec_order<span class="hl opt">)</span> <span class="hl kwa">from</span> dtrios_results_ch.combine<span class="hl opt">(</span> spec_order_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;BBAA_sign_ld05.csv&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;Dmin_sign_ld05.csv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> dtrios_signif_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/dstats.R</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
