---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 10) Analysis VIII (admixture)

## Summary

The population recombination rate is estimated within the [**nextflow**](https://www.nextflow.io/) script `analysis_admixture.nf` (located under `$BASE_DIR/nf/10_analysis_admixture/`), which runs on the XX data set.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:500px; margin:auto;">

</div>

## Details of `analysis_admixture.nf`

### Data preparation

The nextflow script starts by opening the genotype data.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 10.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 10.2</span>
<span class="hl slc">// Set different k values for the admixture analysis</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl num">2</span>..15 <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> k_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 10.3</span>
<span class="hl slc">// load Fst outlier regions</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/plugin/poptrees/outlier.bed&quot;</span><span class="hl opt">)</span>
	.splitCsv<span class="hl opt">(</span>header<span class="hl opt">:</span>true<span class="hl opt">,</span> sep<span class="hl opt">:</span><span class="hl str">&quot;</span><span class="hl esc">\t</span><span class="hl str">&quot;</span><span class="hl opt">)</span>
	.map<span class="hl opt">{</span> row <span class="hl opt">-&gt; [</span> chrom<span class="hl opt">:</span>row.chrom<span class="hl opt">,</span> start<span class="hl opt">:</span>row.start<span class="hl opt">,</span> end<span class="hl opt">:</span>row.end<span class="hl opt">,</span> gid<span class="hl opt">:</span>row.gid <span class="hl opt">] }</span>
	.combine<span class="hl opt">(</span> vcf_ch <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_admx <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 10.4</span>
<span class="hl slc">// subset genotypes to the outlier region and reformat</span>
<span class="hl kwa">process</span> plink12 <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_plink12&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${grouping.gid}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> grouping <span class="hl opt">),</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_admx

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> grouping <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hapmap.*.ped&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hapmap.*.map&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hapmap.*.nosex&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;pop.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> admx_plink

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo -e &quot;CHROM</span><span class="hl esc">\\</span><span class="hl str">tSTART</span><span class="hl esc">\\</span><span class="hl str">tEND&quot; &gt; outl.bed</span>
<span class="hl str">	echo -e &quot;</span><span class="hl ipl">${grouping.chrom}</span><span class="hl str"></span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl ipl">${grouping.start}</span><span class="hl str"></span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl ipl">${grouping.end}</span><span class="hl str">&quot; &gt;&gt; outl.bed</span>
<span class="hl str"></span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep -v &quot;tor</span><span class="hl esc">\\</span><span class="hl str">|tab</span><span class="hl esc">\\</span><span class="hl str">|flo&quot; | \</span>
<span class="hl str">		awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;\$1}&#39; | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">t.*</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)</span><span class="hl esc">\\</span><span class="hl str">(...</span><span class="hl esc">\\</span><span class="hl str">)\$/</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">2/g&#39; &gt; pop.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--keep pop.txt \</span>
<span class="hl str">		--bed outl.bed \</span>
<span class="hl str">		--plink \</span>
<span class="hl str">		--out admx_plink</span>
<span class="hl str"></span>
<span class="hl str">	plink \</span>
<span class="hl str">		--file admx_plink \</span>
<span class="hl str">		--recode12 \</span>
<span class="hl str">		--out hapmap.</span><span class="hl ipl">${grouping.gid}</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 10.5</span>
<span class="hl slc">// combine genoutype subsets with k values</span>
admx_prep  <span class="hl opt">=</span> k_ch.combine<span class="hl opt">(</span> admx_plink <span class="hl opt">)</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 10.6</span>
<span class="hl slc">// run admixture</span>
<span class="hl kwa">process</span> admixture_all <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g4h_admixture_all&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/admixture/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${grouping.gid}</span><span class="hl str">.</span><span class="hl ipl">${k}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> k <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> grouping <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> ped <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> map <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> nosex <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop <span class="hl opt">)</span> <span class="hl kwa">from</span> admx_prep

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;dummy&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.out&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.Q&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> admx_log

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	mv</span> <span class="hl ipl">${pop}</span> <span class="hl str">pop.</span><span class="hl ipl">${grouping.gid}</span><span class="hl str">.</span><span class="hl ipl">${k}</span><span class="hl str">.txt</span>
<span class="hl str">	admixture --cv</span> <span class="hl ipl">${ped} ${k}</span> <span class="hl str">| tee log.</span><span class="hl ipl">${grouping.gid}</span><span class="hl str">.</span><span class="hl ipl">${k}</span><span class="hl str">.out</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
</div>

Finally, we are done with XX.

---
