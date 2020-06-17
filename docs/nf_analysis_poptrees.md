---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 7) Analysis V (poptrees)

## Summary

The population recombination rate is estimated within the [**nextflow**](https://www.nextflow.io/) script `analysis_XX.nf` (located under `$BASE_DIR/nf/0x_analysis_xx/`), which runs on the XX data set.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:500px; margin:auto;">

</div>

## Details of `analysis_xx.nf`

### Data preparation

The nextflow script starts by opening the genotype data.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 7.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_fst <span class="hl opt">}</span>

<span class="hl slc">// git 7.2</span>
<span class="hl slc">// load all possible population pairs</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/plugin/poptrees/all_crosses.tsv&quot;</span><span class="hl opt">)</span>
	.splitCsv<span class="hl opt">(</span>header<span class="hl opt">:</span>true<span class="hl opt">,</span> sep<span class="hl opt">:</span><span class="hl str">&quot;</span><span class="hl esc">\t</span><span class="hl str">&quot;</span><span class="hl opt">)</span>
	.map<span class="hl opt">{</span> row <span class="hl opt">-&gt; [</span> pop1<span class="hl opt">:</span>row.pop1<span class="hl opt">,</span> pop2<span class="hl opt">:</span>row.pop2 <span class="hl opt">] }</span>
	.set<span class="hl opt">{</span> crosses_ch <span class="hl opt">}</span>

<span class="hl slc">// git 7.3</span>
<span class="hl slc">// open the focal Fst outlier regions</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/plugin/poptrees/outlier.bed&quot;</span><span class="hl opt">)</span>
	.splitCsv<span class="hl opt">(</span>header<span class="hl opt">:</span>true<span class="hl opt">,</span> sep<span class="hl opt">:</span><span class="hl str">&quot;</span><span class="hl esc">\t</span><span class="hl str">&quot;</span><span class="hl opt">)</span>
	.map<span class="hl opt">{</span> row <span class="hl opt">-&gt; [</span> chrom<span class="hl opt">:</span>row.chrom<span class="hl opt">,</span> start<span class="hl opt">:</span>row.start<span class="hl opt">,</span> end<span class="hl opt">:</span>row.end<span class="hl opt">,</span> gid<span class="hl opt">:</span>row.gid <span class="hl opt">] }</span>
	.combine<span class="hl opt">(</span> vcf_fst <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> crosses_ch <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> crosses_vcf <span class="hl opt">}</span>

<span class="hl slc">// git 7.4</span>
<span class="hl slc">// compute the average Fst for all possible pair within the outlier region</span>
<span class="hl kwa">process</span> outlier_fst <span class="hl opt">{</span>
		<span class="hl kwb">label</span> <span class="hl str">&quot;L_loc_collect_fst&quot;</span>
		<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/fst/poptree/single&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

		<span class="hl kwb">input</span><span class="hl opt">:</span>
		<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> grouping <span class="hl opt">),</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> cross_pop <span class="hl opt">)</span> <span class="hl kwa">from</span> crosses_vcf

		<span class="hl kwb">output</span><span class="hl opt">:</span>
		<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> grouping.gid <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> cross_pop <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.fst.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> outlier_fst_gid_ch

		<span class="hl kwb">script</span><span class="hl opt">:</span>
		<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">		echo -e &quot;CHROM</span><span class="hl esc">\\</span><span class="hl str">tSTART</span><span class="hl esc">\\</span><span class="hl str">tEND&quot; &gt; outl.bed</span>
<span class="hl str">		echo -e &quot;</span><span class="hl ipl">${grouping.chrom}</span><span class="hl str"></span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl ipl">${grouping.start}</span><span class="hl str"></span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl ipl">${grouping.end}</span><span class="hl str">&quot; &gt;&gt; outl.bed</span>
<span class="hl str"></span>
<span class="hl str">		vcfsamplenames ${vcf[0]} | \</span>
<span class="hl str">			grep &quot;</span><span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">&quot; &gt; pop1.pop</span>
<span class="hl str"></span>
<span class="hl str">			vcfsamplenames ${vcf[0]} | \</span>
<span class="hl str">				grep &quot;</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">&quot; &gt; pop2.pop</span>
<span class="hl str"></span>
<span class="hl str">		vcftools --gzvcf ${vcf[0]} \</span>
<span class="hl str">			--bed outl.bed \</span>
<span class="hl str">			--keep pop1.pop \</span>
<span class="hl str">			--keep pop2.pop \</span>
<span class="hl str">			--weir-fst-pop pop1.pop \</span>
<span class="hl str">			--weir-fst-pop pop2.pop \</span>
<span class="hl str">			--stdout 2&gt;</span> <span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">.50k.log | \</span>
<span class="hl str">			gzip &gt;</span> <span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">.fst.tsv.gz</span>
<span class="hl str"></span>
<span class="hl str">		mFST=\$(grep &quot;Weir and Cockerham mean Fst estimate:&quot;</span> <span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">.50k.log | sed &#39;s/Weir and Cockerham mean Fst estimate: //&#39;)</span>
<span class="hl str">		wFST=\$(grep &quot;Weir and Cockerham weighted Fst estimate:&quot;</span> <span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">.50k.log | sed &#39;s/Weir and Cockerham weighted Fst estimate: //&#39;)</span>
<span class="hl str"></span>
<span class="hl str">		echo -e &quot;</span><span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str"></span><span class="hl esc">\\</span><span class="hl str">t\$mFST</span><span class="hl esc">\\</span><span class="hl str">t\$wFST&quot; &gt;</span> <span class="hl ipl">${cross_pop.pop1}</span><span class="hl str">-</span><span class="hl ipl">${cross_pop.pop2}</span><span class="hl str">.</span><span class="hl ipl">${grouping.gid}</span><span class="hl str">.fst.tsv</span>
<span class="hl str">		&quot;&quot;&quot;</span>
	<span class="hl opt">}</span>

<span class="hl slc">// git 7.5</span>
<span class="hl slc">// collect all population pairs within each region and compile Fst table</span>
<span class="hl kwa">process</span> outlier_fst_collect <span class="hl opt">{</span>
		<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_outlier_fst&quot;</span>
		<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/fst/poptree/summary&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

		<span class="hl kwb">input</span><span class="hl opt">:</span>
		<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> gid <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> cross_pop <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> fst <span class="hl opt">)</span> <span class="hl kwa">from</span> outlier_fst_gid_ch.groupTuple<span class="hl opt">()</span>

		<span class="hl kwb">output</span><span class="hl opt">:</span>
		<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${gid}</span><span class="hl str">.fst.all.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> outlier_fst_collect_ch

		<span class="hl kwb">script</span><span class="hl opt">:</span>
		<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">		echo -e &quot;run</span><span class="hl esc">\\</span><span class="hl str">tmean_fst</span><span class="hl esc">\\</span><span class="hl str">tweighted_fst&quot; &gt;</span> <span class="hl ipl">${gid}</span><span class="hl str">.fst.all.tsv</span>
<span class="hl str">		cat *.fst.tsv &gt;&gt;</span> <span class="hl ipl">${gid}</span><span class="hl str">.fst.all.tsv</span>
<span class="hl str">		&quot;&quot;&quot;</span>
	<span class="hl opt">}</span>

<span class="hl slc">// git 7.6</span>
<span class="hl slc">// neighbour joining happens within the R script R/fig/plot_F4.R</span>
</code>
</pre>
</div>
</div>

Finally, we are done with XX.

---
