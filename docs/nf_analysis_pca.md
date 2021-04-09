---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 73) Analysis V (Principal Component Analysis)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/07_analysis_pca
source ../sh/nextflow_alias.sh
nf_run_pca
```

## Summary

The <span style="color:red;">...</span> are computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_pca.nf` (located under `$BASE_DIR/nf/07_analysis_pca/`).
It takes the <span style="color:red;">...</span> and computes <span style="color:red;">...</span>.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:800px; margin:auto;">

</div>

## Details of `analysis_pca.nf`

### Setup

The nextflow script starts by <span style="color:red;">...</span>

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 7.1</span>
<span class="hl slc">// prepare subset modes (whole genome vs non-diverged regions)</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;whg&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;subset_non_diverged&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> subset_type_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.2</span>
<span class="hl slc">// load table with differentiation outlier regions</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span> <span class="hl str">&quot;../../2_analysis/summaries/fst_outliers_998.tsv&quot;</span> <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> outlier_tab <span class="hl opt">}</span>
</code>
</pre>
</div>


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.3</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> outlier_tab <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> subset_type_ch <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_ch <span class="hl opt">}</span>
</code>
</pre>
</div>


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.4</span>
<span class="hl slc">// depending on subset mode, subset vcf</span>
<span class="hl kwa">process</span> subset_vcf_divergence_based <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_subset_divergence&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_tab <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz.tbi&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_locations<span class="hl opt">,</span> vcf_all_samples_pca <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	if [ &quot;</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot; == &quot;subset_non_diverged&quot; ];then</span>
<span class="hl str">		awk -v OFS=&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; &#39;{print \$2,\$3,\$4}&#39;</span> <span class="hl ipl">${outlier_tab}</span> <span class="hl str">&gt; diverged_regions.bed </span>
<span class="hl str">		SUBSET=&quot;--exclude-bed diverged_regions.bed&quot;</span>
<span class="hl str">	else</span>
<span class="hl str">		SUBSET=&quot;&quot;</span>
<span class="hl str">	fi</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		\$SUBSET \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt;</span> <span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix</span> <span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.5</span>
<span class="hl slc">// prepare location channel for separate pcas</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;bel&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;hon&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;pan&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> locations_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.6</span>
<span class="hl slc">// define location specific sepcies set</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;ind&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;may&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> bel_spec1_ch<span class="hl opt">;</span> bel_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;abe&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;gum&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;ran&quot;</span><span class="hl opt">], [</span><span class="hl num">6</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> hon_spec1_ch<span class="hl opt">;</span> hon_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> pan_spec1_ch<span class="hl opt">;</span> pan_spec2_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.7</span>
<span class="hl slc">// attach genotypes to location channel</span>
locations_ch
	.combine<span class="hl opt">(</span> vcf_locations <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_location_combo <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.8</span>
<span class="hl slc">// subset vcf by location</span>
<span class="hl kwa">process</span> subset_vcf_by_location <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_subset_vcf&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_location_combo

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.pop&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_loc_pca <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| \</span>
<span class="hl str">		grep</span> <span class="hl ipl">${loc}</span> <span class="hl str">| \</span>
<span class="hl str">		grep -v tor | \</span>
<span class="hl str">		grep -v tab &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop</span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--keep</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop \</span>
<span class="hl str">		--mac 3 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | gzip &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// PCA section</span>
<span class="hl slc">// -----------</span>
<span class="hl slc">// git 7.9</span>
<span class="hl slc">// run pca by location</span>
<span class="hl kwa">process</span> pca_location <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g15h_pca_location&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../figures/pca&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> <span class="hl opt">,</span> pattern<span class="hl opt">:</span> <span class="hl str">&quot;*.pdf&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/pca&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> <span class="hl opt">,</span> pattern<span class="hl opt">:</span> <span class="hl str">&quot;*.gz&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_loc_pca

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.prime_pca.pdf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.pca.pdf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.exp_var.txt.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.scores.txt.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> pca_loc_out

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;\$1}&#39;</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">t.*</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)\$/</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">2/g&#39; &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R</span> <span class="hl ipl">${vcf} ${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt 6</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.10</span>
<span class="hl slc">// run pca for global data set</span>
<span class="hl kwa">process</span> pca_all <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g15h_pca_all&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../figures/pca&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> <span class="hl opt">,</span> pattern<span class="hl opt">:</span> <span class="hl str">&quot;*.pdf&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/pca&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> <span class="hl opt">,</span> pattern<span class="hl opt">:</span> <span class="hl str">&quot;*.txt.gz&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../1_genotyping/4_phased/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> <span class="hl opt">,</span> pattern<span class="hl opt">:</span> <span class="hl str">&quot;*.vcf.gz&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_all_samples_pca

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.prime_pca.pdf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.pca.pdf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.exp_var.txt.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.scores.txt.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> pca_all_out
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> vcf_hamlets_only
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> vcf_multi_fst

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# complete PCA, all samples ------------</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| \</span>
<span class="hl str">		awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;\$1}&#39; | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">t.*</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)\$/</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">2/g&#39; &gt; all.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R</span> <span class="hl ipl">${vcf}</span> <span class="hl str">all.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt 6</span>
<span class="hl str"></span>
<span class="hl str">	# PCA without outgroups ---------------</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| \</span>
<span class="hl str">		grep -v &quot;abe</span><span class="hl esc">\\</span><span class="hl str">|gum</span><span class="hl esc">\\</span><span class="hl str">|ind</span><span class="hl esc">\\</span><span class="hl str">|may</span><span class="hl esc">\\</span><span class="hl str">|nig</span><span class="hl esc">\\</span><span class="hl str">|pue</span><span class="hl esc">\\</span><span class="hl str">|ran</span><span class="hl esc">\\</span><span class="hl str">|uni</span><span class="hl esc">\\</span><span class="hl str">|flo&quot; &gt; outgroup.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| \</span>
<span class="hl str">		grep &quot;abe</span><span class="hl esc">\\</span><span class="hl str">|gum</span><span class="hl esc">\\</span><span class="hl str">|ind</span><span class="hl esc">\\</span><span class="hl str">|may</span><span class="hl esc">\\</span><span class="hl str">|nig</span><span class="hl esc">\\</span><span class="hl str">|pue</span><span class="hl esc">\\</span><span class="hl str">|ran</span><span class="hl esc">\\</span><span class="hl str">|uni</span><span class="hl esc">\\</span><span class="hl str">|flo&quot; | \</span>
<span class="hl str">		awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;\$1}&#39; | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">t.*</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)\$/</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">2/g&#39; &gt; hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--remove outgroup.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | gzip &gt; hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/vcf2pca.R hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz hamlets_only.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.pop.txt 6</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
