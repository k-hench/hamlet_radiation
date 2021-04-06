---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 12) Analysis X (<i>F<sub>ST</sub></i> Permutation Test)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/12_analysis_fst_signif
source ../sh/nextflow_alias.sh
nf_run_fstsig
```

## Summary

The <span style="color:red;">...</span> are computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_fst_sign.nf` (located under `$BASE_DIR/nf/12_analysis_fst_signif/`).
It takes the <span style="color:red;">...</span> and computes <span style="color:red;">...</span>.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:800px; margin:auto;">

</div>

## Details of `analysis_fst_sign.nf`

### Setup

The nextflow script starts by <span style="color:red;">...</span>

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 12.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.into<span class="hl opt">{</span> vcf_locations<span class="hl opt">;</span> vcf_adapt <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.2</span>
<span class="hl slc">// prepare location channel</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;bel&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;hon&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;pan&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> locations_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.3</span>
<span class="hl slc">// prepare subset modes (whole genome vs non-diverged regions)</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;whg&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;subset_non_diverged&quot;</span><span class="hl opt">)</span>
	.into<span class="hl opt">{</span> subset_type_ch<span class="hl opt">;</span> subset_type_ch2 <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.4</span>
<span class="hl slc">// load table with differentiation outlier regions</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span> <span class="hl str">&quot;../../2_analysis/summaries/fst_outliers_998.tsv&quot;</span> <span class="hl opt">)</span>
	.into<span class="hl opt">{</span> outlier_tab<span class="hl opt">;</span> outlier_tab2 <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.5</span>
<span class="hl slc">// attach genotypes to location</span>
locations_ch
	.combine<span class="hl opt">(</span> vcf_locations <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> outlier_tab <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_location_combo <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.6</span>
<span class="hl slc">// subset vcf by location</span>
<span class="hl kwa">process</span> subset_vcf_by_location <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_subset_vcf&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_tab <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_location_combo

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${loc}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${loc}</span><span class="hl str">.vcf.gz.tbi&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${loc}</span><span class="hl str">.pop&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_tab <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_loc_pair1<span class="hl opt">,</span> vcf_loc_pair2<span class="hl opt">,</span> vcf_loc_pair3 <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep</span> <span class="hl ipl">${loc}</span> <span class="hl str">| \</span>
<span class="hl str">		grep -v tor | \</span>
<span class="hl str">		grep -v tab &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--keep</span> <span class="hl ipl">${loc}</span><span class="hl str">.pop \</span>
<span class="hl str">		--mac 3 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	</span>
<span class="hl str">	tabix</span> <span class="hl ipl">${loc}</span><span class="hl str">.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.7</span>
<span class="hl slc">// define location specific sepcies set</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;ind&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;may&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> bel_spec1_ch<span class="hl opt">;</span> bel_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;abe&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;gum&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;ran&quot;</span><span class="hl opt">], [</span><span class="hl num">6</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> hon_spec1_ch<span class="hl opt">;</span> hon_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> pan_spec1_ch<span class="hl opt">;</span> pan_spec2_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.8</span>
<span class="hl slc">// prepare pairwise fsts</span>
<span class="hl slc">// ------------------------------</span>
<span class="hl com">/* (create all possible species pairs depending on location</span>
<span class="hl com">   and combine with genotype subset (for the respective location))*/</span>
<span class="hl slc">// ------------------------------</span>
<span class="hl com">/* channel content after joinig:</span>
<span class="hl com">  set [0:val(loc), 1:file(vcf), 2:file( vcfidx ), 3:file(pop), 4:file( outlier_tab ), 5:val(spec1), 6:val(spec2)]*/</span>
<span class="hl slc">// ------------------------------</span>
bel_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;bel&quot;</span> <span class="hl opt">)</span>
	.join<span class="hl opt">(</span> vcf_loc_pair1 <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>bel_spec1_ch<span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>bel_spec2_ch<span class="hl opt">)</span>
	.filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">7</span><span class="hl opt">] }</span>
	.map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">3</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">,</span><span class="hl num">8</span><span class="hl opt">]}</span>
hon_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;hon&quot;</span> <span class="hl opt">)</span>
	.join<span class="hl opt">(</span> vcf_loc_pair2 <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>hon_spec1_ch<span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>hon_spec2_ch<span class="hl opt">)</span>
	.filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">7</span><span class="hl opt">] }</span>
	.map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">3</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">,</span><span class="hl num">8</span><span class="hl opt">]}</span>
pan_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;pan&quot;</span> <span class="hl opt">)</span>
	.join<span class="hl opt">(</span> vcf_loc_pair3 <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>pan_spec1_ch<span class="hl opt">)</span>
	.combine<span class="hl opt">(</span>pan_spec2_ch<span class="hl opt">)</span>
	.filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">7</span><span class="hl opt">] }</span>
	.map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">3</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">,</span><span class="hl num">8</span><span class="hl opt">]}</span>
bel_pairs_ch.concat<span class="hl opt">(</span> hon_pairs_ch<span class="hl opt">,</span> pan_pairs_ch  <span class="hl opt">)</span>.set <span class="hl opt">{</span> all_fst_pairs_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.9</span>
<span class="hl slc">// run fst on actual populations</span>
<span class="hl kwa">process</span> fst_run <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g1h_fst_run&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_tab <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> all_fst_pairs_ch.combine<span class="hl opt">(</span> subset_type_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_random_fst_a00.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> rand_header_ch
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;col1.pop&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;prep.pop&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> rand_body_ch

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
<span class="hl str">		--stdout | bgzip &gt;</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	tabix</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	echo -e &quot;0000</span><span class="hl esc">\t</span><span class="hl str">real_pop&quot; &gt; idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz | \</span>
<span class="hl str">		awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;substr(\$1, length(\$1)-5, length(\$1))}&#39;  &gt; prep.pop</span>
<span class="hl str">	grep</span> <span class="hl ipl">${spec1} ${pop}</span> <span class="hl str">&gt; pop1.txt</span>
<span class="hl str">	grep</span> <span class="hl ipl">${spec2} ${pop}</span> <span class="hl str">&gt; pop2.txt</span>
<span class="hl str">	</span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${loc}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz \</span>
<span class="hl str">		--weir-fst-pop pop1.txt \</span>
<span class="hl str">		--weir-fst-pop pop2.txt \</span>
<span class="hl str">		--stdout 2&gt; fst.log 1&gt; tmp.txt</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;^Weir&quot; fst.log | sed &#39;s/.* //&#39; | paste - - &gt; fst.tsv</span>
<span class="hl str">	echo -e &quot;idx</span><span class="hl esc">\\</span><span class="hl str">ttype</span><span class="hl esc">\\</span><span class="hl str">tmean_fst</span><span class="hl esc">\\</span><span class="hl str">tweighted_fst&quot; &gt;</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">_random_fst_a00.tsv</span>
<span class="hl str">	paste idx.txt fst.tsv &gt;&gt;</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">_random_fst_a00.tsv</span>
<span class="hl str"></span>
<span class="hl str">	rm fst.tsv fst.log pop1.txt pop2.txt tmp.txt idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{print \$1}&#39; prep.pop &gt; col1.pop</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.10</span>
<span class="hl slc">// create indexes for permutation itteration</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;0&#39;</span>..<span class="hl str">&#39;9&#39;</span><span class="hl opt">))</span>
	.map<span class="hl opt">{</span> <span class="hl str">&quot;0&quot;</span> <span class="hl opt">+</span> it <span class="hl opt">}</span> <span class="hl slc">//.into{ sub_pre_ch; sub_pre_ch2 }</span>
	.into<span class="hl opt">{</span> singles_ch<span class="hl opt">;</span> tens_ch <span class="hl opt">}</span>

singles_ch
	.combine<span class="hl opt">(</span>tens_ch<span class="hl opt">)</span>
	.map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">]+</span>it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">] }</span>
	.toSortedList<span class="hl opt">()</span>
	.flatten<span class="hl opt">()</span>
	.into<span class="hl opt">{</span> sub_pre_ch<span class="hl opt">;</span> sub_pre_ch2 <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.11</span>
<span class="hl slc">// for each itteration run fst on 100</span>
<span class="hl slc">// permutations of population assignment</span>
<span class="hl kwa">process</span> random_bodies <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g6h_fst_run&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> col1 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> prepop <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> pre <span class="hl opt">)</span> <span class="hl kwa">from</span> rand_body_ch.combine<span class="hl opt">(</span>sub_pre_ch<span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span><span class="hl str">&quot;*_random_fst_b</span><span class="hl ipl">${pre}</span><span class="hl str">.tsv&quot;</span><span class="hl opt">)</span> <span class="hl kwa">into</span> rand_body_out_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	for k in {00..99}; do</span>
<span class="hl str">	echo &quot;Iteration_&quot;\$k</span>
<span class="hl str">	echo -e &quot;</span><span class="hl ipl">${pre}</span><span class="hl str">\$k</span><span class="hl esc">\t</span><span class="hl str">random&quot; &gt; idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{print \$2}&#39;</span> <span class="hl ipl">${prepop}</span> <span class="hl str">| shuf &gt; col2.pop # premutation happens here</span>
<span class="hl str">	paste</span> <span class="hl ipl">${col1}</span> <span class="hl str">col2.pop &gt; rand.pop</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">\$&quot; rand.pop &gt; r_pop1.pop</span>
<span class="hl str">	grep &quot;</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">\$&quot; rand.pop &gt; r_pop2.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--weir-fst-pop r_pop1.pop \</span>
<span class="hl str">		--weir-fst-pop r_pop2.pop \</span>
<span class="hl str">		--stdout  2&gt; fst.log 1&gt; tmp.txt</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;^Weir&quot; fst.log | sed &#39;s/.* //&#39; | paste - - &gt; fst.tsv</span>
<span class="hl str">	paste idx.txt fst.tsv &gt;&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst_b</span><span class="hl ipl">${pre}</span><span class="hl str">.tsv</span>
<span class="hl str"></span>
<span class="hl str">	rm fst.tsv fst.log rand.pop col2.pop r_pop1.pop r_pop2.pop tmp.txt </span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.12</span>
<span class="hl slc">// collect all itterations and compile</span>
<span class="hl slc">// output for each population pair</span>
<span class="hl kwa">process</span> compile_random_results <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_compile_rand&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/fst_signif/random&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> body <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> head <span class="hl opt">)</span> <span class="hl kwa">from</span> rand_body_out_ch.groupTuple<span class="hl opt">()</span>.join<span class="hl opt">(</span>rand_header_ch<span class="hl opt">,</span> remainder<span class="hl opt">:</span> true<span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span><span class="hl str">&quot;</span><span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv.gz&quot;</span><span class="hl opt">)</span> <span class="hl kwa">into</span> random_lists_result

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	cat</span> <span class="hl ipl">${head}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	cat</span> <span class="hl ipl">${body}</span> <span class="hl str">&gt;&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	gzip</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// -----------------------------------------</span>
<span class="hl slc">// repeat the same procedure for adaptation</span>
<span class="hl slc">// (permuting location within species)</span>

<span class="hl slc">// git 12.13</span>
<span class="hl slc">// prepare species channel</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> species_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.14</span>
<span class="hl slc">// define location set</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;bel&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;hon&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;pan&quot;</span><span class="hl opt">]])</span>.into<span class="hl opt">{</span> locations_ch_1<span class="hl opt">;</span>locations_ch_2 <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.15</span>
<span class="hl slc">// create location pairs</span>
locations_ch_1
	.combine<span class="hl opt">(</span>locations_ch_2<span class="hl opt">)</span>
	.filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">2</span><span class="hl opt">] }</span>
	.map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">3</span><span class="hl opt">]}</span>
	.combine<span class="hl opt">(</span> species_ch <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> vcf_adapt <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> outlier_tab2 <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> subset_type_ch2 <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_location_combo_adapt <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.16</span>
<span class="hl slc">// collapsed analog to git 12.6 &amp; 9</span>
<span class="hl slc">// subset vcf by species and</span>
<span class="hl slc">// run fst on actual populations</span>
<span class="hl kwa">process</span> fst_run_adapt <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g1h_fst_run&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc2 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> vcf_indx<span class="hl opt">) ,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_tab <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> subset_type <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_location_combo_adapt

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}${loc1}</span><span class="hl str">-</span><span class="hl ipl">${spec}${loc2}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_random_fst_a00.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> rand_header_adapt_ch
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}${loc1}</span><span class="hl str">-</span><span class="hl ipl">${spec}${loc2}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;col1.pop&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;prep.pop&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> rand_body_adapt_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep</span> <span class="hl ipl">${spec}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${spec}</span><span class="hl str">.pop</span>
<span class="hl str"></span>
<span class="hl str">	if [ &quot;</span><span class="hl ipl">${subset_type}</span><span class="hl str">&quot; == &quot;subset_non_diverged&quot; ];then</span>
<span class="hl str">		awk -v OFS=&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; &#39;{print \$2,\$3,\$4}&#39;</span> <span class="hl ipl">${outlier_tab}</span> <span class="hl str">&gt; diverged_regions.bed </span>
<span class="hl str">		SUBSET=&quot;--exclude-bed diverged_regions.bed&quot;</span>
<span class="hl str">	else</span>
<span class="hl str">		SUBSET=&quot;&quot;</span>
<span class="hl str">	fi</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		\$SUBSET \</span>
<span class="hl str">		--keep</span> <span class="hl ipl">${spec}</span><span class="hl str">.pop \</span>
<span class="hl str">		--mac 3 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt;</span> <span class="hl ipl">${spec}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	tabix</span> <span class="hl ipl">${spec}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	echo -e &quot;0000</span><span class="hl esc">\t</span><span class="hl str">real_pop&quot; &gt; idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${spec}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz | \</span>
<span class="hl str">		awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;substr(\$1, length(\$1)-5, length(\$1))}&#39;  &gt; prep.pop</span>
<span class="hl str">	grep</span> <span class="hl ipl">${loc1} ${spec}</span><span class="hl str">.pop &gt; pop1.txt</span>
<span class="hl str">	grep</span> <span class="hl ipl">${loc2} ${spec}</span><span class="hl str">.pop &gt; pop2.txt</span>
<span class="hl str">	</span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${spec}</span><span class="hl str">.</span><span class="hl ipl">${subset_type}</span><span class="hl str">.vcf.gz \</span>
<span class="hl str">		--weir-fst-pop pop1.txt \</span>
<span class="hl str">		--weir-fst-pop pop2.txt \</span>
<span class="hl str">		--stdout 2&gt; fst.log 1&gt; tmp.txt</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;^Weir&quot; fst.log | sed &#39;s/.* //&#39; | paste - - &gt; fst.tsv</span>
<span class="hl str">	echo -e &quot;idx</span><span class="hl esc">\\</span><span class="hl str">ttype</span><span class="hl esc">\\</span><span class="hl str">tmean_fst</span><span class="hl esc">\\</span><span class="hl str">tweighted_fst&quot; &gt;</span> <span class="hl ipl">${spec}${loc1}</span><span class="hl str">-</span><span class="hl ipl">${spec}${loc2}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">_random_fst_a00.tsv</span>
<span class="hl str">	paste idx.txt fst.tsv &gt;&gt;</span> <span class="hl ipl">${spec}${loc1}</span><span class="hl str">-</span><span class="hl ipl">${spec}${loc2}</span><span class="hl str">_</span><span class="hl ipl">${subset_type}</span><span class="hl str">_random_fst_a00.tsv</span>
<span class="hl str"></span>
<span class="hl str">	rm fst.tsv fst.log pop1.txt pop2.txt tmp.txt idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{print \$1}&#39; prep.pop &gt; col1.pop</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.17</span>
<span class="hl slc">// for each itteration run fst on 100</span>
<span class="hl slc">// permutations of location assignment</span>
<span class="hl kwa">process</span> random_bodies_adapt <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g6h_fst_run&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> col1 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> prepop <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> pre <span class="hl opt">)</span> <span class="hl kwa">from</span> rand_body_adapt_ch.combine<span class="hl opt">(</span>sub_pre_ch2<span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span><span class="hl str">&quot;*_random_fst_b</span><span class="hl ipl">${pre}</span><span class="hl str">.tsv&quot;</span><span class="hl opt">)</span> <span class="hl kwa">into</span> rand_body_out_adapt_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	for k in {00..99}; do</span>
<span class="hl str">	echo &quot;Iteration_&quot;\$k</span>
<span class="hl str">	echo -e &quot;</span><span class="hl ipl">${pre}</span><span class="hl str">\$k</span><span class="hl esc">\t</span><span class="hl str">random&quot; &gt; idx.txt</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{print \$2}&#39;</span> <span class="hl ipl">${prepop}</span> <span class="hl str">| shuf &gt; col2.pop # premutation happens here</span>
<span class="hl str">	paste</span> <span class="hl ipl">${col1}</span> <span class="hl str">col2.pop &gt; rand.pop</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;</span><span class="hl ipl">${spec}${loc1}</span><span class="hl str">\$&quot; rand.pop &gt; r_pop1.pop</span>
<span class="hl str">	grep &quot;</span><span class="hl ipl">${spec}${loc2}</span><span class="hl str">\$&quot; rand.pop &gt; r_pop2.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--weir-fst-pop r_pop1.pop \</span>
<span class="hl str">		--weir-fst-pop r_pop2.pop \</span>
<span class="hl str">		--stdout  2&gt; fst.log 1&gt; tmp.txt</span>
<span class="hl str"></span>
<span class="hl str">	grep &quot;^Weir&quot; fst.log | sed &#39;s/.* //&#39; | paste - - &gt; fst.tsv</span>
<span class="hl str">	paste idx.txt fst.tsv &gt;&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst_b</span><span class="hl ipl">${pre}</span><span class="hl str">.tsv</span>
<span class="hl str"></span>
<span class="hl str">	rm fst.tsv fst.log rand.pop col2.pop r_pop1.pop r_pop2.pop tmp.txt </span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 12.18</span>
<span class="hl slc">// collect all itterations and compile</span>
<span class="hl slc">// output for each location pair</span>
<span class="hl kwa">process</span> compile_random_results_adapt <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_compile_rand&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/fst_signif/random/adapt&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> run <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> body <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> head <span class="hl opt">)</span> <span class="hl kwa">from</span> rand_body_out_adapt_ch.groupTuple<span class="hl opt">()</span>.join<span class="hl opt">(</span>rand_header_adapt_ch<span class="hl opt">,</span> remainder<span class="hl opt">:</span> true<span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span><span class="hl str">&quot;</span><span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv.gz&quot;</span><span class="hl opt">)</span> <span class="hl kwa">into</span> random_lists_adapt_result

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	cat</span> <span class="hl ipl">${head}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	cat</span> <span class="hl ipl">${body}</span> <span class="hl str">&gt;&gt;</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	gzip</span> <span class="hl ipl">${run}</span><span class="hl str">_random_fst.tsv</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
