---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 14) Analysis XII (revPoMo Phylogenies)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/14_analysis_revpomo
source ../sh/nextflow_alias.sh
nf_run_revpomo
```

## Summary

The <span style="color:red;">...</span> are computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_revpomo.nf` (located under `$BASE_DIR/nf/14_analysis_revpomo/`).
It takes the <span style="color:red;">...</span> and computes <span style="color:red;">...</span>.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:800px; margin:auto;">

</div>

## Details of `analysis_revpomo.nf`

### Setup

The nextflow script starts by <span style="color:red;">...</span>

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 14.1</span>
<span class="hl slc">// Open the SNP data set</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_snps_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.2</span>
<span class="hl slc">// Prepare LG channel for vcf subsetting</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;01&#39;</span>..<span class="hl str">&#39;09&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;10&#39;</span>..<span class="hl str">&#39;19&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;20&#39;</span>..<span class="hl str">&#39;24&#39;</span><span class="hl opt">))</span>
	.map<span class="hl opt">{</span> <span class="hl str">&quot;LG&quot;</span> <span class="hl opt">+</span> it <span class="hl opt">}</span>
	.set<span class="hl opt">{</span> lg_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.3</span>
<span class="hl slc">// Subset snp data set by LG</span>
<span class="hl kwa">process</span> subset_snps_by_lg <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_subset_lg&quot;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span> <span class="hl opt">(</span> lg <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_snps_ch.map<span class="hl opt">{ [</span>it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">]</span>.minus<span class="hl opt">(</span><span class="hl str">&quot;.vcf&quot;</span><span class="hl opt">),</span> it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">]]}</span>.combine<span class="hl opt">(</span> lg_ch <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz*&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> vcf_snps_lg_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--chr</span> <span class="hl ipl">${lg}</span> <span class="hl str">\</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix</span> <span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.4</span>
<span class="hl slc">// Open the allBP data set (will be expanded x 24 LGs)</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG*.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_allbp_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.5</span>
<span class="hl slc">// Open the pre-defined window positions</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/windows_1kb.bed.gz&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> windows_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.6</span>
<span class="hl slc">// Subset ALL vcf files (also allBP) by missingnes (max. 10%)</span>
<span class="hl kwa">process</span> filter_vcf_missingnes <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g6h_filter_vcf&quot;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_snps_lg_ch.concat<span class="hl opt">(</span> vcf_allbp_ch <span class="hl opt">)</span>.map<span class="hl opt">{ [</span>it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">]</span>.minus<span class="hl opt">(</span><span class="hl str">&quot;.vcf&quot;</span><span class="hl opt">),</span> it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">]]}</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> vcfId <span class="hl opt">),</span>  <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_single_ind.vcf.gz*&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> vcf_snps_filterd_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--max-missing 0.9 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		sed &quot;s/&lt;NON_REF&gt;/./g&quot; | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_filtered.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_filtered.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_filtered.vcf.gz | head -n 1 &gt; first_ind.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_filtered.vcf.gz \</span>
<span class="hl str">		--keep first_ind.txt \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_single_ind.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_single_ind.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.7</span>
<span class="hl slc">// Coverage of SNPs vcf for SNPdensity, allBP for Ns</span>
<span class="hl kwa">process</span> compute_coverage <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_50g3h_coverage&quot;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/revPoMo/coverage&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> window <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_snps_filterd_ch.combine<span class="hl opt">(</span> windows_ch <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">_cov.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> coverage_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	LG=\$( echo</span> <span class="hl ipl">${vcfId}</span> <span class="hl str">| sed &#39;s/.*</span><span class="hl esc">\\</span><span class="hl str">(LG[0-9]</span><span class="hl esc">\\</span><span class="hl str">{2</span><span class="hl esc">\\</span><span class="hl str">}</span><span class="hl esc">\\</span><span class="hl str">)/</span><span class="hl esc">\\</span><span class="hl str">1/&#39; )</span>
<span class="hl str">	echo -e &quot;CHROM</span><span class="hl esc">\\</span><span class="hl str">tSTART</span><span class="hl esc">\\</span><span class="hl str">tEND&quot; &gt; windows_1kb.\$LG.bed</span>
<span class="hl str"></span>
<span class="hl str">	zcat</span> <span class="hl ipl">${window}</span> <span class="hl str">| \</span>
<span class="hl str">		grep \$LG &gt;&gt; windows_1kb.\$LG.bed</span>
<span class="hl str">	</span>
<span class="hl str">	gzip windows_1kb.\$LG.bed</span>
<span class="hl str"></span>
<span class="hl str">	bedtools coverage \</span>
<span class="hl str">		-a windows_1kb.\$LG.bed.gz \</span>
<span class="hl str">		-b</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		-counts  &gt;</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_cov.tsv</span>
<span class="hl str">	</span>
<span class="hl str">	gzip</span> <span class="hl ipl">${vcfId}</span><span class="hl str">_cov.tsv</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.8</span>
<span class="hl slc">// Compile summary table</span>
<span class="hl kwa">process</span> compile_window_stats <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_window_stats&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/revPoMo/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> windows <span class="hl opt">)</span> <span class="hl kwa">from</span> coverage_ch.collect<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;window_stats.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> final_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	#!/usr/bin/env Rscript</span>
<span class="hl str"></span>
<span class="hl str">	library(tidyverse)</span>
<span class="hl str"></span>
<span class="hl str">	data_SNPs &lt;- 1:24 %&gt;% </span>
<span class="hl str">					str_pad(width = 2, pad = 0) %&gt;%</span>
<span class="hl str">					str_c(&quot;phased_mac2.LG&quot;, ., &quot;_cov.tsv.gz&quot;) %&gt;%</span>
<span class="hl str">					map_dfr(.f = function(file){</span>
<span class="hl str">							read_tsv(file, </span>
<span class="hl str">								col_names = c(&quot;CHROM&quot;, &quot;START&quot;, &quot;END&quot;, &quot;COV_SNP&quot;)) %&gt;%</span>
<span class="hl str">							filter(COV_SNP &gt; 0 )</span>
<span class="hl str">								} ) </span>
<span class="hl str">					</span>
<span class="hl str">	data_allBPs &lt;- 1:24 %&gt;% </span>
<span class="hl str">					str_pad(width = 2, pad = 0) %&gt;%</span>
<span class="hl str">					str_c(&quot;filterd.allBP.LG&quot;, ., &quot;_cov.tsv.gz&quot;) %&gt;%</span>
<span class="hl str">					map_dfr(.f = function(file){</span>
<span class="hl str">							read_tsv(file, </span>
<span class="hl str">								col_names = c(&quot;CHROM&quot;, &quot;START&quot;, &quot;END&quot;, &quot;COV_ALL&quot;)) %&gt;%</span>
<span class="hl str">							filter(COV_ALL &gt; 0 )</span>
<span class="hl str">							}</span>
<span class="hl str">						)</span>
<span class="hl str"></span>
<span class="hl str">	data &lt;- data_SNPs %&gt;%</span>
<span class="hl str">		left_join(data_allBPs, by = c(CHROM = &quot;CHROM&quot;, START = &quot;START&quot;, END = &quot;END&quot;)) %&gt;%</span>
<span class="hl str">		mutate(SNP_density = COV_SNP/ COV_ALL, </span>
<span class="hl str">				REL_COV =  COV_ALL/ (END-START))</span>
<span class="hl str">	</span>
<span class="hl str">	data %&gt;% write_tsv(&quot;window_stats.tsv.gz&quot;)</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
