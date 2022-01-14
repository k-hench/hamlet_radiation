---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---



# (git 13) Analysis XI (Whole Genome Phylogenies)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/13_analysis_phylo_whg
nextflow run analysis_phylo_whg.nf
```

## Summary

The whole genome phylogenies can be reconstructed within the [**nextflow**](https://www.nextflow.io/) script `analysis_phylo_whg.nf` (located under `$BASE_DIR/nf/analysis_phylo_whg/`).

## Details of `analysis_phylo_whg.nf`

> This part of the analysis was actually manged manually and not via `nextflow`. 
> We still report the analysis as a `.nf` script as we believe this is a cleaner and more concise report of the conducted analysis.

### Setup

The nextflow script starts by opening the genotype data and feeding it into two different streams.

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow

<span class="hl slc">// git 13.1</span>
<span class="hl slc">// Open the SNP data set</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.into<span class="hl opt">{</span> vcf_snps_ch<span class="hl opt">;</span> vcf_snps_ch2 <span class="hl opt">}</span>
</code>
</pre>
</div>

We also open the *allBP* data set (genotypes of SNPS + invariant sites).


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.2</span>
<span class="hl slc">// Open the allBP data set (will be expanded x 24 LGs)</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/filterd.allBP.non_ref.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_allbp_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Then, we initialize the different window sizes considered (1-50 kb).


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.3</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( [</span> <span class="hl num">1</span><span class="hl opt">,</span> <span class="hl num">5</span><span class="hl opt">,</span> <span class="hl num">10</span><span class="hl opt">,</span> <span class="hl num">50</span> <span class="hl opt">] )</span>
	.set<span class="hl opt">{</span> window_size_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Then we create one bed file per windowsize, containing the window positions.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.4</span>
<span class="hl slc">// Compile summary table</span>
<span class="hl kwa">process</span> segment_windows <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_loc_slice_windows&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/window_stats/windows/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">)</span> <span class="hl kwa">from</span> window_size_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;windows_</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb.bed.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> windows_ch<span class="hl opt">,</span> windows_ch2<span class="hl opt">,</span> windows_ch3 <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	#!/usr/bin/env Rscript</span>
<span class="hl str">	library(hypogen)</span>
<span class="hl str"></span>
<span class="hl str">	x_bp &lt;-</span> <span class="hl ipl">${kb_size}</span> <span class="hl str">* 1000</span>
<span class="hl str">	</span>
<span class="hl str">	window_non_overlap &lt;- function(CHROM, LENGTH, windowsize = x_bp ){</span>
<span class="hl str">	  tibble(CHROM = CHROM, </span>
<span class="hl str">	         START = seq(from = 1, to = LENGTH, by = windowsize) - 1, # (produces overlap of one which is needed for bedtools)</span>
<span class="hl str">	         END = lead(START, default = LENGTH) ) }</span>
<span class="hl str"></span>
<span class="hl str">	hypo_karyotype %&gt;%</span>
<span class="hl str">	  select(CHROM, LENGTH) %&gt;%</span>
<span class="hl str">	  pmap_dfr(window_non_overlap) %&gt;%</span>
<span class="hl str">	  write_tsv(file = &quot;windows_</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb.bed.gz&quot;)</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

We drop the Serranid samples (the outgroup) from the genotypes.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.5</span>
<span class="hl kwa">process</span> filter_hamlets <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_outgroup_drop&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> vcfId <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_snps_ch2

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;hamlets_only&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hamlets_only.vcf.gz*&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> vcf_hamlets_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo -e &quot;20478tabhon</span><span class="hl esc">\\</span><span class="hl str">n28393torpan</span><span class="hl esc">\\</span><span class="hl str">ns_tort_3torpan&quot; &gt; outgroup.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools  --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--remove outgroup.pop \</span>
<span class="hl str">		--mac 1 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		bgzip &gt; hamlets_only.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix hamlets_only.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then, we compute the coverage (n of SNPs) within each window.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.6</span>
<span class="hl slc">// Coverage of SNPs vcf for SNPdensity, allBP for Ns</span>
<span class="hl kwa">process</span> compute_coverage <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_140g1h_coverage&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/window_stats/coverages/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl slc">//set vcfId, file( vcf ), val( kb_size ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )</span>
	<span class="hl kwa">set</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> vcfId <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> window <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_snps_ch.concat<span class="hl opt">(</span> vcf_hamlets_ch <span class="hl opt">)</span>.map<span class="hl opt">{ [</span>it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">]</span>.minus<span class="hl opt">(</span><span class="hl str">&quot;.vcf&quot;</span><span class="hl opt">),</span> it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">]]}</span>.combine<span class="hl opt">(</span> windows_ch <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> coverage_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	bedtools coverage \</span>
<span class="hl str">		-a</span> <span class="hl ipl">${window}</span> <span class="hl str">\</span>
<span class="hl str">		-b</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		-counts | gzip &gt;</span> <span class="hl ipl">${vcfId}</span><span class="hl str">.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

To split the genotypes by linkage group, we initialize a LG-channel...


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.7</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;01&#39;</span>..<span class="hl str">&#39;09&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;10&#39;</span>..<span class="hl str">&#39;19&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;20&#39;</span>..<span class="hl str">&#39;24&#39;</span><span class="hl opt">) )</span>
	.set<span class="hl opt">{</span> lg_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

...and split the genotypes.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.8</span>
<span class="hl kwa">process</span> subset_allBP <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_140g10h_coverage&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/non_ref_byLG/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl slc">//set vcfId, file( vcf ), val( kb_size ), file( window ) from vcf_snps_filterd_ch.combine( windows_ch )</span>
	<span class="hl kwa">set</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">)</span> <span class="hl kwa">from</span> lg_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;filterd.allBP.non_ref.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;filterd.allBP.non_ref.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz.tbi&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> allbp_non_ref_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered//filterd.allBP.non_ref.vcf.gz \</span>
<span class="hl str">		--chr LG</span><span class="hl ipl">${lg}</span> <span class="hl str">\</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt; filterd.allBP.non_ref.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	tabix filterd.allBP.non_ref.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then we also compute the coverage (in terms of callable sites, not just SNPs) within windows.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.9</span>
<span class="hl kwa">process</span> compute_coverage_allBP <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_140g1h_coverage&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/window_stats/coverages/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> tbi <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> window <span class="hl opt">)</span> <span class="hl kwa">from</span> allbp_non_ref_ch.combine<span class="hl opt">(</span> windows_ch2 <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;filterd.allBP.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> coverage_allbp_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo -e &quot;CHROM</span><span class="hl esc">\\</span><span class="hl str">tSTART</span><span class="hl esc">\\</span><span class="hl str">tEND&quot; &gt; bed_LG</span><span class="hl ipl">${lg}</span><span class="hl str">.bed</span>
<span class="hl str"></span>
<span class="hl str">	zgrep &quot;LG</span><span class="hl ipl">${lg}</span><span class="hl str">&quot;</span> <span class="hl ipl">${window}</span> <span class="hl str">&gt;&gt; bed_LG</span><span class="hl ipl">${lg}</span><span class="hl str">.bed</span>
<span class="hl str">	gzip bed_LG</span><span class="hl ipl">${lg}</span><span class="hl str">.bed</span>
<span class="hl str"></span>
<span class="hl str">	bedtools coverage \</span>
<span class="hl str">		-a bed_LG</span><span class="hl ipl">${lg}</span><span class="hl str">.bed.gz \</span>
<span class="hl str">		-b</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		-counts | gzip &gt; filterd.allBP.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

So, now we can merge the coverage information for each window.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.10</span>
<span class="hl slc">// Compile summary table</span>
<span class="hl kwa">process</span> complie_window_stats <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_windows_stats&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/window_stats/window_stats/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> windows <span class="hl opt">)</span> <span class="hl kwa">from</span> coverage_ch.concat<span class="hl opt">(</span> coverage_allbp_ch <span class="hl opt">)</span>.groupTuple<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;window_stats.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> window_out_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	#!/usr/bin/env Rscript</span>
<span class="hl str"></span>
<span class="hl str">	library(tidyverse)</span>
<span class="hl str"></span>
<span class="hl str">	data_SNPs &lt;- read_tsv(&quot;phased_mac2.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz&quot;,</span>
<span class="hl str">						  col_names = c(&quot;CHROM&quot;, &quot;START&quot;, &quot;END&quot;, &quot;COV_SNP&quot;))</span>
<span class="hl str"></span>
<span class="hl str">	data_HYP &lt;- read_tsv(&quot;hamlets_only.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_cov.tsv.gz&quot;,</span>
<span class="hl str">						  col_names = c(&quot;CHROM&quot;, &quot;START&quot;, &quot;END&quot;, &quot;COV_HYP&quot;))</span>
<span class="hl str"></span>
<span class="hl str">	all_bp_files &lt;- dir(pattern = &quot;filterd.allBP.LG*&quot;)</span>
<span class="hl str"></span>
<span class="hl str">	data_allBPs &lt;- map_dfr(all_bp_files, .f = function(x){read_tsv(x, col_names = c(&quot;CHROM&quot;, &quot;START&quot;, &quot;END&quot;, &quot;COV_ALL&quot;))})</span>
<span class="hl str"></span>
<span class="hl str">	data &lt;- data_SNPs %&gt;%</span>
<span class="hl str">		left_join(data_HYP, by = c(CHROM = &quot;CHROM&quot;, START = &quot;START&quot;, END = &quot;END&quot;))  %&gt;%</span>
<span class="hl str">		left_join(data_allBPs, by = c(CHROM = &quot;CHROM&quot;, START = &quot;START&quot;, END = &quot;END&quot;))  %&gt;%</span>
<span class="hl str">		filter(COV_ALL &gt; 0 ) %&gt;%</span>
<span class="hl str">		mutate(SNP_density = round(COV_SNP/ COV_ALL, 2), </span>
<span class="hl str">		REL_COV =  round(COV_ALL/ (END-START), 2))</span>
<span class="hl str">	</span>
<span class="hl str">	write_tsv(x = data, file = &quot;window_stats.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb.tsv.gz&quot;)</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

From this point onward, the analysis was actually done manually (we still format it according to nextflow for clarity).


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// form here on, this pipeline was not actually run using</span>
<span class="hl slc">// nextflow, but managed manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>
</code>
</pre>
</div>

Now, the subset of windows to be considered for the phylogeny was sampled.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.11</span>
<span class="hl slc">// Subset to 5000 random windows</span>
<span class="hl kwa">process</span> subset_windows <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_loc_subset_windows&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> windows <span class="hl opt">)</span> <span class="hl kwa">from</span> window_out_ch.filter<span class="hl opt">({</span> it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">] ==</span> <span class="hl num">5</span> <span class="hl opt">})</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb_size <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;5000x_</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_v1.bed&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> window_subset_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	#!/usr/bin/env Rscript</span>
<span class="hl str">	library(hypogen)   # attaches tidyverse automatically</span>
<span class="hl str"></span>
<span class="hl str">	### Draw windows with coverage (sequence contiguity) and SNP-based cutoffs</span>
<span class="hl str"></span>
<span class="hl str">	stats5 &lt;- read_tsv(&quot;window_stats.</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb.tsv.gz&quot;)</span>
<span class="hl str"></span>
<span class="hl str">	set.seed(64)</span>
<span class="hl str"></span>
<span class="hl str">	selected_windows &lt;- stats5 %&gt;%</span>
<span class="hl str">	filter(REL_COV &gt;= 0.80 &amp;</span>
<span class="hl str">			COV_HYP &gt;= 50) %&gt;%</span>
<span class="hl str">	sample_n(5000) %&gt;%</span>
<span class="hl str">	arrange(CHROM, START)</span>
<span class="hl str"></span>
<span class="hl str">	write_tsv(selected_windows %&gt;% select(CHROM, START, END), file = &quot;5000x_</span><span class="hl ipl">${kb_size}</span><span class="hl str">kb_v1.bed&quot;)</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

To parallelize the analysis, several looping-channels are opened.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.12</span>
<span class="hl slc">// index for sub-channels</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl num">1</span>..50 <span class="hl opt">)</span>.into<span class="hl opt">{</span> loop50_idx_ch1<span class="hl opt">;</span> loop50_idx_ch2<span class="hl opt">;</span> loop50_idx_ch3<span class="hl opt">;</span> loop50_idx_ch4 <span class="hl opt">}</span>
</code>
</pre>
</div>

For the selected windows, the genotypes are extracted.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.13</span>
<span class="hl slc">// Extract windows from all BP</span>
<span class="hl kwa">process</span> extract_windows <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_extract_windows&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> idx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> windows <span class="hl opt">)</span> <span class="hl kwa">from</span> loop50_idx_ch1.combine<span class="hl opt">(</span> window_out_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_v1_all.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> window_extract_loop
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	PER_TASK=100</span>
<span class="hl str"></span>
<span class="hl str">	START_NUM=\$(( (</span><span class="hl ipl">${idx}</span> <span class="hl str">- 1) * \$PER_TASK + 1 ))</span>
<span class="hl str">	END_NUM=\$((</span> <span class="hl ipl">${idx}</span> <span class="hl str">* \$PER_TASK ))</span>
<span class="hl str"></span>
<span class="hl str">	echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, which will do runs \$START_NUM to \$END_NUM</span>
<span class="hl str"></span>
<span class="hl str">	for (( run=\$START_NUM ; run&lt;=END_NUM ; run++ ))</span>
<span class="hl str">		do</span>
<span class="hl str">		echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, run number \$run</span>
<span class="hl str">		chr=\$(awk -v line=&quot;\$run&quot; &#39;BEGIN { FS = &quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; } ; NR==line+1 { print \$1 }&#39; 5000x_5kb_v1.bed)</span>
<span class="hl str">		sta=\$(awk -v line=&quot;\$run&quot; &#39;BEGIN { FS = &quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; } ; NR==line+1 { print \$2 }&#39; 5000x_5kb_v1.bed)</span>
<span class="hl str">		end=\$(awk -v line=&quot;\$run&quot; &#39;BEGIN { FS = &quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; } ; NR==line+1 { print \$3 }&#39; 5000x_5kb_v1.bed)</span>
<span class="hl str">		printf -v i &quot;%04d&quot; \$run</span>
<span class="hl str"></span>
<span class="hl str">		vcftools \</span>
<span class="hl str">			--gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/byLG/filterd.allBP.&quot;\$chr&quot;.vcf.gz \</span>
<span class="hl str">			--chr &quot;\$chr&quot; \</span>
<span class="hl str">			--from-bp &quot;\$sta&quot; \</span>
<span class="hl str">			--to-bp &quot;\$end&quot; \</span>
<span class="hl str">			--recode --stdout | \</span>
<span class="hl str">			grep -v &quot;##&quot; | bgzip &gt; window_&quot;\$i&quot;_v1_all.vcf.gz</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Again, a channel is created for parallelizing...


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.14</span>
<span class="hl slc">// index for sub-channels</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl num">1</span>..20 <span class="hl opt">)</span>.set<span class="hl opt">{</span> loop20_idx_ch <span class="hl opt">}</span>
<span class="hl slc">// bundle the distributed sub-channels</span>
window_extract_loop.collect<span class="hl opt">()</span>.map<span class="hl opt">{ [</span> it <span class="hl opt">] }</span>.set<span class="hl opt">{</span> window_extract_ch1<span class="hl opt">,</span> window_extract_ch2 <span class="hl opt">}</span>
</code>
</pre>
</div>

...and again, the outgroup samples are removed from the genoptypes.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.15</span>
<span class="hl slc">// Remove samples (all / noS data sets)</span>
<span class="hl kwa">process</span> remove_samples <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_remove_samples&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> idx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> loop20_idx_ch.combine<span class="hl opt">(</span> window_extract_ch1 <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_noS.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> samples_removed_loop
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	PER_TASK=250</span>
<span class="hl str"></span>
<span class="hl str">	START_NUM=\$(( (</span><span class="hl ipl">${idx}</span> <span class="hl str">- 1) * \$PER_TASK + 1 ))</span>
<span class="hl str">	END_NUM=\$((</span> <span class="hl ipl">${idx}</span> <span class="hl str">* \$PER_TASK ))</span>
<span class="hl str"></span>
<span class="hl str">	echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, which will do runs \$START_NUM to \$END_NUM</span>
<span class="hl str"></span>
<span class="hl str">	for (( run=\$START_NUM ; run&lt;=END_NUM ; run++ ))</span>
<span class="hl str">		do</span>
<span class="hl str"></span>
<span class="hl str">		echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, run number \$run</span>
<span class="hl str">		</span>
<span class="hl str">		printf -v i &quot;%04d&quot; \$run</span>
<span class="hl str"></span>
<span class="hl str">		vcftools \</span>
<span class="hl str">			--gzvcf window_&quot;\$i&quot;*.vcf.gz \</span>
<span class="hl str">			--remove samples_Serr.txt \</span>
<span class="hl str">			--recode \</span>
<span class="hl str">			--stdout | bgzip &gt; window_&quot;\$i&quot;_v1_noS.vcf.gz</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
<span class="hl slc">// bundle the distributed sub-channels</span>
samples_removed_loop.collect<span class="hl opt">()</span>.map<span class="hl opt">{ [</span> it <span class="hl opt">] }</span>.set<span class="hl opt">{</span> samples_removed_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Then, we translate the genotypes from vcf to fasta format.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.16</span>
<span class="hl slc">// Convert to Fasta (IUPAC-encoded)</span>
<span class="hl kwa">process</span> fasta_convert <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_convert_to_fasta&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> idx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf_noS <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf_all <span class="hl opt">)</span> <span class="hl kwa">from</span> loop50_idx_ch2.combine<span class="hl opt">(</span> samples_removed_ch <span class="hl opt">)</span>.combine<span class="hl opt">(</span> window_extract_ch2 <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> fasta_convert_loop
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	PER_TASK=100</span>
<span class="hl str"></span>
<span class="hl str">	START_NUM=\$(( (</span><span class="hl ipl">${idx}</span> <span class="hl str">- 1) * \$PER_TASK + 1 ))</span>
<span class="hl str">	END_NUM=\$((</span> <span class="hl ipl">${idx}</span> <span class="hl str">* \$PER_TASK ))</span>
<span class="hl str"></span>
<span class="hl str">	echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, which will do runs \$START_NUM to \$END_NUM</span>
<span class="hl str"></span>
<span class="hl str">	for (( run=\$START_NUM ; run&lt;=END_NUM ; run++ ))</span>
<span class="hl str">		do</span>
<span class="hl str">		echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, run number \$run</span>
<span class="hl str">		printf -v i &quot;%04d&quot; \$run</span>
<span class="hl str">		for j in all noS</span>
<span class="hl str">			do</span>
<span class="hl str">				bgzip -cd window_&quot;\$i&quot;_v1_&quot;\$j&quot;.vcf.gz | vcf-to-tab &gt; window_&quot;\$i&quot;_v1_&quot;\$j&quot;.tab</span>
<span class="hl str"></span>
<span class="hl str">				perl \$SFTWR/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i window_&quot;\$i&quot;_v1_&quot;\$j&quot;.tab &gt; window_&quot;\$i&quot;_v1_&quot;\$j&quot;.fas</span>
<span class="hl str"></span>
<span class="hl str">				rm window_&quot;\$i&quot;_v1_&quot;\$j&quot;.tab*</span>
<span class="hl str">		done</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
<span class="hl slc">// bundle the distributed sub-channels</span>
fasta_convert_loop.collect<span class="hl opt">()</span>.map<span class="hl opt">{ [</span> it <span class="hl opt">] }</span>.set<span class="hl opt">{</span> fasta_convert_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

The fasta sequences are aligned.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.17</span>
<span class="hl slc">// Align sequences in windows</span>
<span class="hl kwa">process</span> fasta_align <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_fasta_align&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> idx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> fa <span class="hl opt">)</span> <span class="hl kwa">from</span> loop50_idx_ch3.combine<span class="hl opt">(</span> fasta_convert_ch <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.aln&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> fasta_align_loop
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	PER_TASK=100</span>
<span class="hl str"></span>
<span class="hl str">	START_NUM=\$(( (</span><span class="hl ipl">${idx}</span> <span class="hl str">- 1) * \$PER_TASK + 1 ))</span>
<span class="hl str">	END_NUM=\$((</span> <span class="hl ipl">${idx}</span> <span class="hl str">* \$PER_TASK ))</span>
<span class="hl str"></span>
<span class="hl str">	echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, which will do runs \$START_NUM to \$END_NUM</span>
<span class="hl str"></span>
<span class="hl str">	for (( run=</span><span class="hl ipl">${idx}</span> <span class="hl str">; run&lt;=END_NUM ; run++ ))</span>
<span class="hl str">		do</span>
<span class="hl str">		echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, run number \$run</span>
<span class="hl str">		printf -v i &quot;%04d&quot; \$run</span>
<span class="hl str">		for j in all noS</span>
<span class="hl str">			do</span>
<span class="hl str">				mafft --auto window_&quot;\$i&quot;_v1_&quot;\$j&quot;.fas &gt; &quot;window_&quot;\$i&quot;_v1_&quot;\$j&quot;.aln</span>
<span class="hl str">		done</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
<span class="hl slc">// bundle the distributed sub-channels</span>
fasta_align_loop.collect<span class="hl opt">()</span>.map<span class="hl opt">{ [</span> it <span class="hl opt">] }</span>.set<span class="hl opt">{</span> fasta_align_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Next, the trees are inferred within each window.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.18</span>
<span class="hl slc">// Infer local trees</span>
<span class="hl kwa">process</span> local_trees <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_local_trees&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> idx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> aln <span class="hl opt">)</span> <span class="hl kwa">from</span> loop50_idx_ch4.combine<span class="hl opt">(</span> fasta_align_ch <span class="hl opt">)</span>
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.treefile&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> local_trees_loop
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	PER_TASK=100</span>
<span class="hl str"></span>
<span class="hl str">	START_NUM=$(( (</span><span class="hl ipl">${idx}</span><span class="hl str">- 1) * \$PER_TASK + 1 ))</span>
<span class="hl str">	END_NUM=$((</span> <span class="hl ipl">${idx}</span> <span class="hl str">* \$PER_TASK ))</span>
<span class="hl str"></span>
<span class="hl str">	echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, which will do runs \$START_NUM to \$END_NUM</span>
<span class="hl str"></span>
<span class="hl str">	for (( run=</span><span class="hl ipl">${idx}</span> <span class="hl str">; run&lt;=END_NUM ; run++ ))</span>
<span class="hl str">		do</span>
<span class="hl str">		echo This is task</span> <span class="hl ipl">${idx}</span><span class="hl str">, run number \$run</span>
<span class="hl str">		printf -v i &quot;%04d&quot; \$run</span>
<span class="hl str">		for j in all noS</span>
<span class="hl str">			do</span>
<span class="hl str">				iqtree2 -s window_&quot;\$i&quot;_v1_&quot;\$j&quot;.aln  --prefix locus_&quot;\$i&quot;_v1_&quot;\$j&quot; -o PL17_160floflo -T 1</span>
<span class="hl str">		done</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
<span class="hl slc">// bundle the distributed sub-channels</span>
local_trees_loop.collect<span class="hl opt">()</span>.map<span class="hl opt">{ [</span> it <span class="hl opt">] }</span>.set<span class="hl opt">{</span> local_trees_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Finally, the individual trees are summarized.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.19</span>
<span class="hl slc">// Calculate summary tree</span>
<span class="hl kwa">process</span> summary_trees <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_summary_trees&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/astral/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> tree <span class="hl opt">)</span> <span class="hl kwa">from</span> local_trees_ch
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;astral*.tre&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;astral*.log&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> summary_trees_ch
	
	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	cat ./*_noS.treefile &gt; genetrees_5000x_5kb_v1_noS.tre</span>
<span class="hl str">	java -jar \$SFTWR/ASTRAL_5.7.5/astral.5.7.5.jar -i genetrees_5000x_5kb_v1_noS.tre -o astral_5000x_5kb_v1_noS.tre 2&gt; astral_5000x_5kb_v1_noS.log</span>
<span class="hl str"></span>
<span class="hl str">	cat ./*_all.treefile &gt; genetrees_5000x_5kb_v1_all.tre</span>
<span class="hl str">	java -jar \$SFTWR/ASTRAL_5.7.5/astral.5.7.5.jar -i genetrees_5000x_5kb_v1_all.tre -o astral_5000x_5kb_v1_all.tre 2&gt; astral_5000x_5kb_v1_all.log</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
