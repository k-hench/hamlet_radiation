---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 16) Analysis XIV (Identity by Descent)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/16_analysis_ibd
nextflow run analysis_ibd.nf -c ../../nextflow.config -resume
```

## Summary

Identity by descent is computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_ibd.nf` (located under `$BASE_DIR/nf/16_analysis_ibd/`).
It takes the phased genotypes and computes the IBD segments.
Below is an overview of the steps involved in the analysis.

## Details of `analysis_ibd.nf`

### Setup

The nextflow script starts by opening the phased genotypes.

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// This pipeline includes the analysis calculating the pair-wise IBD</span>

<span class="hl slc">// git 16.1</span>
<span class="hl slc">// load genotypes</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> genotypes_raw_ch <span class="hl opt">}</span>
</code>
</pre>
</div>


Then, outgroup samples are removed.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 16.2</span>
<span class="hl slc">// drop outgroups</span>
<span class="hl kwa">process</span> drop_outgroups <span class="hl opt">{</span>
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> genotypes_raw_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;phased_mac2.no_outgroup.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> genotypes_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep &quot;tor</span><span class="hl esc">\\</span><span class="hl str">|tab</span><span class="hl esc">\\</span><span class="hl str">|flo&quot; &gt; outrgr.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--remove outrgr.pop \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt; phased_mac2.no_outgroup.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

The IDB analysis is run for several IBD fragment size thresholds, so here we are initializing the different thresholds.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 16.3</span>
<span class="hl slc">// Set IBD fragment sizes</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">([[</span> <span class="hl num">25000</span><span class="hl opt">,</span> <span class="hl num">10000</span><span class="hl opt">,</span> <span class="hl num">7</span> <span class="hl opt">],</span>
	       <span class="hl opt">[</span> <span class="hl num">15000</span><span class="hl opt">,</span> <span class="hl num">7500</span><span class="hl opt">,</span> <span class="hl num">10</span> <span class="hl opt">],</span>
	       <span class="hl opt">[</span> <span class="hl num">10000</span><span class="hl opt">,</span> <span class="hl num">5000</span><span class="hl opt">,</span> <span class="hl num">8</span> <span class="hl opt">]])</span>
	.set<span class="hl opt">{</span> seq_sizes_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Also, the IBD analysis is additionally run with specific parts of the genome excluded - here the different exclusion models are initialized.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 16.4</span>
<span class="hl slc">// Set filter mode</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">([[</span><span class="hl str">&quot;direct&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;&quot;</span><span class="hl opt">],</span>
	       <span class="hl opt">[</span><span class="hl str">&quot;bed&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;95&quot;</span><span class="hl opt">]])</span>
	.set<span class="hl opt">{</span> filtermode_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Truffle is run to compute the IBD segments.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 16.5</span>
<span class="hl slc">// run truffle</span>
<span class="hl kwa">process</span> run_truffle <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/ibd/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sz1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sz2 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sz3 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> mode <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> excluding <span class="hl opt">)</span> <span class="hl kwa">from</span> genotypes_ch.combine<span class="hl opt">(</span> seq_sizes_ch <span class="hl opt">)</span>.combine<span class="hl opt">(</span> filtermode_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sz3 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd.tsv&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> truffle_result

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	if<span class="hl opt">(</span> mode <span class="hl opt">==</span> <span class="hl str">&#39;direct&#39;</span> <span class="hl opt">)</span>
		<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">		truffle \</span>
<span class="hl str">			--vcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">			--segments \</span>
<span class="hl str">			--nofiltering \</span>
<span class="hl str">			--ibs1markers</span> <span class="hl ipl">${sz1}</span> <span class="hl str">\</span>
<span class="hl str">			--ibs2markers</span> <span class="hl ipl">${sz2}</span> <span class="hl str">\</span>
<span class="hl str">			--out no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span> <span class="hl str">\</span>
<span class="hl str">			--cpu 8</span>
<span class="hl str">		</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd.tsv</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments.tsv</span>
<span class="hl str">		&quot;&quot;&quot;</span>
	else if<span class="hl opt">(</span> mode <span class="hl opt">==</span> <span class="hl str">&#39;filter&#39;</span> <span class="hl opt">)</span>
		<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">		vcftools --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">			--not-chr</span> <span class="hl ipl">${excluding}</span> <span class="hl str">\</span>
<span class="hl str">			--recode \</span>
<span class="hl str">			--stdout | bgzip &gt; tmp.vcf.gz </span>
<span class="hl str">		</span>
<span class="hl str">		truffle \</span>
<span class="hl str">			--vcf tmp.vcf.gz  \</span>
<span class="hl str">			--segments \</span>
<span class="hl str">			--nofiltering \</span>
<span class="hl str">			--ibs1markers</span> <span class="hl ipl">${sz1}</span> <span class="hl str">\</span>
<span class="hl str">			--ibs2markers</span> <span class="hl ipl">${sz2}</span> <span class="hl str">\</span>
<span class="hl str">			--out no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span> <span class="hl str">\</span>
<span class="hl str">			--cpu 8</span>
<span class="hl str">		</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd.tsv</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments.tsv</span>
<span class="hl str">		rm tmp.vcf.gz </span>
<span class="hl str">		&quot;&quot;&quot;</span>
	else if<span class="hl opt">(</span> mode <span class="hl opt">==</span> <span class="hl str">&#39;bed&#39;</span> <span class="hl opt">)</span>
		<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">		vcftools --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">			--exclude-bed ../../../../../ressources/plugin/idb_above_</span><span class="hl ipl">${excluding}</span><span class="hl str">.bed \</span>
<span class="hl str">			--recode \</span>
<span class="hl str">			--stdout | bgzip &gt; tmp.vcf.gz </span>
<span class="hl str">		</span>
<span class="hl str">		truffle \</span>
<span class="hl str">			--vcf tmp.vcf.gz  \</span>
<span class="hl str">			--segments \</span>
<span class="hl str">			--nofiltering \</span>
<span class="hl str">			--ibs1markers</span> <span class="hl ipl">${sz1}</span> <span class="hl str">\</span>
<span class="hl str">			--ibs2markers</span> <span class="hl ipl">${sz2}</span> <span class="hl str">\</span>
<span class="hl str">			--out no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span> <span class="hl str">\</span>
<span class="hl str">			--cpu 8</span>
<span class="hl str">		</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.ibd.tsv</span>
<span class="hl str">		sed &#39;s/^</span><span class="hl esc">\\</span><span class="hl str">s*//g; s/</span><span class="hl esc">\\</span><span class="hl str">s</span><span class="hl esc">\\</span><span class="hl str">+/</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments &gt; no_outgr_</span><span class="hl ipl">${mode}${excluding}</span><span class="hl str">_</span><span class="hl ipl">${sz3}</span><span class="hl str">.segments.tsv</span>
<span class="hl str">		rm tmp.vcf.gz </span>
<span class="hl str">		&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Within R, the genomic coordinates of the IBD segments are converted to cM positions.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 16.6</span>
<span class="hl slc">// convert IBD segments to cM</span>
<span class="hl kwa">process</span> convert_to_cM <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/ibd/cM_converted&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sz3 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> truffle_summary <span class="hl opt">) ,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> truffle_segments <span class="hl opt">)</span> <span class="hl kwa">from</span> truffle_result

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.converted.tsv&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.conv_summary.tsv&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.conv_filterd.tsv&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> cM_result

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	#!/usr/bin/env Rscript</span>
<span class="hl str">	base_dir &lt;- Sys.getenv(&quot;BASE_DIR&quot;)</span>
<span class="hl str"></span>
<span class="hl str">	args &lt;- c(&quot;ressources/recombination/&quot;,</span>
<span class="hl str">          &quot;MAP1cm.txt&quot;, &quot;MAP1bp.txt&quot;,</span>
<span class="hl str">          &quot;MAP2cm.txt&quot;, &quot;MAP2bp.txt&quot;,</span>
<span class="hl str">          &quot;</span><span class="hl ipl">${truffle_segments}</span><span class="hl str">&quot;,</span>
<span class="hl str">          &quot;</span><span class="hl ipl">${truffle_summary}</span><span class="hl str">&quot;)</span>
<span class="hl str">	</span>
<span class="hl str">	renv::activate(base_dir)</span>
<span class="hl str"></span>
<span class="hl str">	library(GenomicOriginsScripts)</span>
<span class="hl str">	library(hypogen)</span>
<span class="hl str">	library(patchwork)</span>
<span class="hl str">	library(plyranges)</span>
<span class="hl str"></span>
<span class="hl str">	rec_path &lt;- str_c(base_dir, as.character(args[1]))</span>
<span class="hl str">	hypo_map1_cm &lt;- as.character(args[2])</span>
<span class="hl str">	hypo_map1_bp &lt;- as.character(args[3])</span>
<span class="hl str">	hypo_map2_cm &lt;- as.character(args[4])</span>
<span class="hl str">	hypo_map2_bp &lt;- as.character(args[5])</span>
<span class="hl str">	segment_file &lt;- as.character(args[6])</span>
<span class="hl str">	summary_file &lt;- as.character(args[7])</span>
<span class="hl str"></span>
<span class="hl str">	truffle_conv &lt;- segment_file %&gt;% str_replace(pattern = &quot;.segments.tsv&quot;, replacement = &quot;.converted.tsv&quot;) %&gt;% str_remove(&quot;.*/&quot;)</span>
<span class="hl str">	truffle_sum &lt;- segment_file %&gt;% str_replace(pattern = &quot;.segments.tsv&quot;, replacement = &quot;.conv_summary.tsv&quot;) %&gt;% str_remove(&quot;.*/&quot;)</span>
<span class="hl str">	truffle_filt &lt;- segment_file %&gt;% str_replace(pattern = &quot;.segments.tsv&quot;, replacement = &quot;.conv_filterd.tsv&quot;) %&gt;% str_remove(&quot;.*/&quot;)</span>
<span class="hl str"></span>
<span class="hl str">	read_maps &lt;- function(cm_file, bp_file){</span>
<span class="hl str">	read_tsv(cm_file) %&gt;% </span>
<span class="hl str">		group_by(LG) %&gt;% </span>
<span class="hl str">		mutate(LGnm = as.roman(LG) %&gt;% as.numeric(),</span>
<span class="hl str">			CHROM = str_c(&quot;LG&quot;, str_pad(LGnm, width = 2, pad = 0)),</span>
<span class="hl str">			cM = ifelse(LG == &quot;VIII&quot;, max(cM)-cM,cM)) %&gt;% </span>
<span class="hl str">		ungroup() %&gt;% </span>
<span class="hl str">		dplyr::select(Loci, cM, CHROM)  %&gt;% </span>
<span class="hl str">		full_join(read_tsv(bp_file) %&gt;% </span>
<span class="hl str">					filter(!(duplicated(Loci) | duplicated(Loci, fromlast = TRUE))),</span>
<span class="hl str">				by = c(Loci = &quot;Loci&quot;, CHROM = &quot;LG&quot;)) %&gt;% </span>
<span class="hl str">		left_join(hypo_chrom_start) %&gt;% </span>
<span class="hl str">		mutate(GPOS = GSTART + bp) %&gt;% </span>
<span class="hl str">		filter(!is.na(GPOS),</span>
<span class="hl str">			!is.na(cM)) %&gt;% </span>
<span class="hl str">		arrange(GPOS)</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	make_lg_seg &lt;- function(lg = &quot;LG08&quot;, n = 31, gmap = gmap1){</span>
<span class="hl str">	data_pos &lt;- tibble(CHROM = rep(lg, n),</span>
<span class="hl str">						start = seq(from = hypo_karyotype\$GSTART[hypo_karyotype\$CHROM == lg],</span>
<span class="hl str">									to = hypo_karyotype\$GEND[hypo_karyotype\$CHROM == lg],</span>
<span class="hl str">									length = n) %&gt;%</span>
<span class="hl str">						floor(),</span>
<span class="hl str">						GSTART = hypo_karyotype\$GSTART[hypo_karyotype\$CHROM == lg]) %&gt;% </span>
<span class="hl str">		mutate(GPOS = start,</span>
<span class="hl str">			start = start - GSTART,</span>
<span class="hl str">			end = start) %&gt;% </span>
<span class="hl str">		as_iranges()</span>
<span class="hl str">	</span>
<span class="hl str">	map_pos &lt;- gmap %&gt;%</span>
<span class="hl str">		filter(CHROM == lg) %&gt;% </span>
<span class="hl str">		attach_end(LG = lg) %&gt;% </span>
<span class="hl str">		mutate(start = lag(bp, default = 0),</span>
<span class="hl str">			start_cM = lag(cM, default = 0)) %&gt;% </span>
<span class="hl str">		dplyr::select(CHROM, start, end = bp, start_cM, end_cM = cM) %&gt;%</span>
<span class="hl str">		mutate(start_bp = start, end_bp = end) %&gt;% </span>
<span class="hl str">		as_iranges()</span>
<span class="hl str">	</span>
<span class="hl str">	list(data = data_pos, map = map_pos)</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	attach_end &lt;- function(data, LG = &quot;LG01&quot;){</span>
<span class="hl str">	data %&gt;% </span>
<span class="hl str">		bind_rows(., </span>
<span class="hl str">				data %&gt;%</span>
<span class="hl str">					filter(row_number() == last(row_number())) %&gt;% </span>
<span class="hl str">					mutate(GPOS = hypo_karyotype\$GEND[hypo_karyotype\$CHROM == CHROM],</span>
<span class="hl str">						bp = hypo_karyotype\$LENGTH[hypo_karyotype\$CHROM == CHROM],</span>
<span class="hl str">						Loci = as.numeric(</span>
<span class="hl str">							str_c(&quot;-99&quot;,</span>
<span class="hl str">								str_remove(string = CHROM, &quot;LG&quot;))</span>
<span class="hl str">						)</span>
<span class="hl str">					))</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	bin_rescaler &lt;- function(bp, start_cM, end_cM, start_bp, end_bp,...){</span>
<span class="hl str">	scales::rescale(x = bp,</span>
<span class="hl str">					to = c(start_cM, end_cM),</span>
<span class="hl str">					from = c(start_bp, end_bp))</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	interpol_data &lt;- function(lg, ...){</span>
<span class="hl str">	data_pair &lt;- make_lg_seg(lg = lg, ...)</span>
<span class="hl str">	</span>
<span class="hl str">	plyranges::join_overlap_inner(data_pair\$data,</span>
<span class="hl str">									data_pair\$map) %&gt;% </span>
<span class="hl str">		as.data.frame() %&gt;% </span>
<span class="hl str">		as_tibble() %&gt;% </span>
<span class="hl str">		dplyr::select(CHROM = CHROM.x, bp = start, GSTART, GPOS, start_cM:end_bp) %&gt;% </span>
<span class="hl str">		mutate(interpol_cM = pmap_dbl(cur_data(), bin_rescaler)) </span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	na_to_zero &lt;- function(x){</span>
<span class="hl str">	x_type &lt;- typeof(x)</span>
<span class="hl str">	if_else(is.na(x), as(0,Class = x_type), x) %&gt;% </span>
<span class="hl str">		as.double() %&gt;% as(Class = x_type)</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	convert_bp_to_cm &lt;- function(data, lg = &quot;LG08&quot;, gmap = gmap1){</span>
<span class="hl str">	gmap_in &lt;- deparse(substitute(gmap))</span>
<span class="hl str">	</span>
<span class="hl str">	data_pos &lt;- data %&gt;% </span>
<span class="hl str">		filter( CHROM == lg ) %&gt;% </span>
<span class="hl str">		dplyr::select(PAIR, TYPE, CHROM, START, END, NMARKERS) %&gt;% </span>
<span class="hl str">		mutate(seg_id = str_c(PAIR,&quot;_&quot;,CHROM,&quot;_&quot;,START)) %&gt;% </span>
<span class="hl str">		pivot_longer(cols = START:END, names_to = &quot;PART&quot;, values_to = &quot;start&quot;) %&gt;%</span>
<span class="hl str">		mutate(end = start) %&gt;% </span>
<span class="hl str">		as_iranges()</span>
<span class="hl str">	</span>
<span class="hl str">	map_pos &lt;- gmap %&gt;%</span>
<span class="hl str">		filter(CHROM == lg) %&gt;%</span>
<span class="hl str">		attach_end(LG = lg) %&gt;%</span>
<span class="hl str">		mutate(start = lag(bp, default = 0) + 1, # avoid overlapping segments - causes duplications in joining</span>
<span class="hl str">			start_cM = lag(cM, default = 0)) %&gt;%</span>
<span class="hl str">		dplyr::select(start, end = bp, start_cM, end_cM = cM) %&gt;%</span>
<span class="hl str">		mutate(start_bp = start, end_bp = end) %&gt;%</span>
<span class="hl str">		as_iranges()</span>
<span class="hl str">	</span>
<span class="hl str">	map_nr &lt;- str_replace(gmap_in, pattern = &quot;gmap&quot;, replacement = &quot;_m&quot;)</span>
<span class="hl str">	</span>
<span class="hl str">	plyranges::join_overlap_inner(data_pos,</span>
<span class="hl str">									map_pos) %&gt;%</span>
<span class="hl str">		as.data.frame() %&gt;%</span>
<span class="hl str">		as_tibble() %&gt;%</span>
<span class="hl str">		dplyr::select(CHROM = CHROM, bp = start, PAIR, TYPE, PART, start_cM:end_bp, seg_id, NMARKERS) %&gt;%</span>
<span class="hl str">		mutate(interpol_cM = pmap_dbl(cur_data(), bin_rescaler)) %&gt;%</span>
<span class="hl str">		pivot_wider(id_cols = c(CHROM,PAIR,TYPE,seg_id, PART,NMARKERS),</span>
<span class="hl str">					values_from = c(bp,interpol_cM),</span>
<span class="hl str">					names_from = PART) %&gt;%</span>
<span class="hl str">		dplyr::select(-seg_id) %&gt;%</span>
<span class="hl str">		mutate(length_bp = bp_END - bp_START,</span>
<span class="hl str">			length_cM = interpol_cM_END - interpol_cM_START) %&gt;% </span>
<span class="hl str">		set_names(value = c(&quot;CHROM&quot;, &quot;PAIR&quot;, &quot;TYPE&quot;, &quot;NMARKERS&quot;, &quot;bp_START&quot;, &quot;bp_END&quot;,</span>
<span class="hl str">							str_c(c(&quot;interpol_cM_START&quot;, &quot;interpol_cM_END&quot;), map_nr),</span>
<span class="hl str">							&quot;length_bp&quot;, str_c(&quot;length_cM&quot;, map_nr)))</span>
<span class="hl str">	}</span>
<span class="hl str"></span>
<span class="hl str">	# actual script -------------------</span>
<span class="hl str">	gmap1 &lt;- read_maps(cm_file = str_c(rec_path, hypo_map1_cm),</span>
<span class="hl str">					bp_file = str_c(rec_path, hypo_map1_bp))</span>
<span class="hl str">	gmap2 &lt;- read_maps(cm_file = str_c(rec_path, hypo_map2_cm),</span>
<span class="hl str">					bp_file = str_c(rec_path, hypo_map2_bp))</span>
<span class="hl str"></span>
<span class="hl str"></span>
<span class="hl str">	lgs &lt;- 1:24 %&gt;%</span>
<span class="hl str">	str_pad(width = 2, pad = 0) %&gt;%</span>
<span class="hl str">	str_c(&quot;LG&quot;,.)</span>
<span class="hl str"></span>
<span class="hl str">	segments_individual_interpol_map1 &lt;- lgs %&gt;% map_dfr(interpol_data, n = 51)</span>
<span class="hl str">	segments_individual_interpol_map2 &lt;- lgs %&gt;% map_dfr(interpol_data, n = 51, gmap = gmap2)</span>
<span class="hl str"></span>
<span class="hl str">	segments_individual &lt;- vroom::vroom(segment_file) %&gt;%</span>
<span class="hl str">	mutate(LENGTH = LENGTH * 10^6,</span>
<span class="hl str">			START = POS * 10^6,</span>
<span class="hl str">			END = START + LENGTH,</span>
<span class="hl str">			PAIR = str_c(ID1, &quot;-&quot;, ID2))</span>
<span class="hl str"></span>
<span class="hl str">	segments_summary &lt;- vroom::vroom(summary_file) %&gt;%</span>
<span class="hl str">	mutate(PAIR = str_c(ID1, &quot;-&quot;, ID2))</span>
<span class="hl str"></span>
<span class="hl str">	bounds_gmap1 &lt;- gmap1 %&gt;% </span>
<span class="hl str">	group_by(CHROM) %&gt;% </span>
<span class="hl str">	filter(cM == max(cM)) %&gt;% </span>
<span class="hl str">	filter(bp == max(bp)) %&gt;% </span>
<span class="hl str">	ungroup() %&gt;% </span>
<span class="hl str">	dplyr::select(CHROM, cM) %&gt;%</span>
<span class="hl str">	mutate(GSTART_cM = cumsum(lag(cM, default = 0)),</span>
<span class="hl str">			GEND_cM = cM + GSTART_cM,</span>
<span class="hl str">			GMID_cM = (GSTART_cM + GEND_cM) / 2,</span>
<span class="hl str">			grp = c(&quot;even&quot;, &quot;odd&quot;)[ 1+row_number() %% 2 ])</span>
<span class="hl str"></span>
<span class="hl str">	bounds_gmap2 &lt;- gmap2 %&gt;% </span>
<span class="hl str">	group_by(CHROM) %&gt;% </span>
<span class="hl str">	filter(cM == max(cM)) %&gt;% </span>
<span class="hl str">	filter(bp == max(bp)) %&gt;% </span>
<span class="hl str">	ungroup() %&gt;% </span>
<span class="hl str">	dplyr::select(CHROM, cM) %&gt;%</span>
<span class="hl str">	mutate(GSTART_cM = cumsum(lag(cM, default = 0)),</span>
<span class="hl str">			GEND_cM = cM + GSTART_cM,</span>
<span class="hl str">			GMID_cM = (GSTART_cM + GEND_cM) / 2,</span>
<span class="hl str">			grp = c(&quot;even&quot;, &quot;odd&quot;)[ 1+row_number() %% 2 ])</span>
<span class="hl str"></span>
<span class="hl str">	hypo_all_starts &lt;- hypo_karyotype %&gt;% </span>
<span class="hl str">	dplyr::select(CHROM, GSTART, GEND) %&gt;% </span>
<span class="hl str">	left_join(bounds_gmap1 %&gt;% </span>
<span class="hl str">				dplyr::select(CHROM, GSTART_cM_m1 = GSTART_cM, GEND_cM_m1 = GEND_cM)) %&gt;% </span>
<span class="hl str">	left_join(bounds_gmap2 %&gt;% </span>
<span class="hl str">				dplyr::select(CHROM, GSTART_cM_m2 = GSTART_cM, GEND_cM_m2 = GEND_cM))</span>
<span class="hl str"></span>
<span class="hl str">	converted_segments &lt;- lgs %&gt;% </span>
<span class="hl str">	map_dfr(convert_bp_to_cm, data = segments_individual) %&gt;% </span>
<span class="hl str">	left_join( lgs %&gt;% </span>
<span class="hl str">				map_dfr(convert_bp_to_cm, data = segments_individual, gmap = gmap2) ) %&gt;%</span>
<span class="hl str">	left_join(hypo_all_starts) %&gt;% </span>
<span class="hl str">	mutate(G_SEG_START = GSTART + bp_START,</span>
<span class="hl str">			G_SEG_END = GSTART + bp_END,</span>
<span class="hl str">			G_SEG_START_cM_m1 = GSTART_cM_m1 + interpol_cM_START_m1,</span>
<span class="hl str">			G_SEG_END_cM_m1 = GSTART_cM_m1 + interpol_cM_END_m1,</span>
<span class="hl str">			G_SEG_START_cM_m2 = GSTART_cM_m2 + interpol_cM_START_m2,</span>
<span class="hl str">			G_SEG_END_cM_m2 = GSTART_cM_m2 + interpol_cM_END_m2)</span>
<span class="hl str"></span>
<span class="hl str">	hypo_cM_length_map1 &lt;- max(hypo_all_starts\$GEND_cM_m1)</span>
<span class="hl str">	hypo_cM_length_map2 &lt;- max(hypo_all_starts\$GEND_cM_m2)</span>
<span class="hl str">	hypo_bp_length &lt;- hypo_karyotype\$GEND[hypo_karyotype\$CHROM == &quot;LG24&quot;]</span>
<span class="hl str"></span>
<span class="hl str">	control &lt;- converted_segments %&gt;%</span>
<span class="hl str">	ungroup() %&gt;%</span>
<span class="hl str">	group_by(PAIR, TYPE) %&gt;%</span>
<span class="hl str">	summarise(seq_length = sum(length_bp),</span>
<span class="hl str">				n_mark = sum(NMARKERS),</span>
<span class="hl str">				cm_length_m1 = sum(length_cM_m1),</span>
<span class="hl str">				cm_length_m2 = sum(length_cM_m2)) %&gt;%</span>
<span class="hl str">	ungroup() %&gt;%</span>
<span class="hl str">	pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:cm_length_m2, values_fill = 0) %&gt;% </span>
<span class="hl str">	left_join(segments_summary, .,  ) %&gt;% </span>
<span class="hl str">	mutate(across(.cols = seq_length_IBD1:cm_length_m2_IBD2, .fns = na_to_zero)) %&gt;% </span>
<span class="hl str">	mutate(IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,</span>
<span class="hl str">			IBD1_manual = n_mark_IBD1 / NMARK,</span>
<span class="hl str">			IBD2_manual = n_mark_IBD2 / NMARK,</span>
<span class="hl str">			icheck_0 = IBD0_manual - IBD0,</span>
<span class="hl str">			icheck_1 = IBD1_manual - IBD1,</span>
<span class="hl str">			icheck_2 = IBD2 - IBD2,</span>
<span class="hl str">			# compile ibd by sequence map</span>
<span class="hl str">			ibd0_bp = (hypo_bp_length - (seq_length_IBD1 + seq_length_IBD2)) / hypo_bp_length,</span>
<span class="hl str">			ibd1_bp = seq_length_IBD1 / hypo_bp_length,</span>
<span class="hl str">			ibd2_bp = seq_length_IBD2 / hypo_bp_length,</span>
<span class="hl str">			# compile ibd by genetic map 1</span>
<span class="hl str">			ibd0_cM_m1 = (hypo_cM_length_map1 - (cm_length_m1_IBD1 + cm_length_m1_IBD2)) / hypo_cM_length_map1,</span>
<span class="hl str">			ibd1_cM_m1 = cm_length_m1_IBD1 / hypo_cM_length_map1,</span>
<span class="hl str">			ibd2_cM_m1 = cm_length_m1_IBD2 / hypo_cM_length_map1,</span>
<span class="hl str">			# compile ibd by genetic map 2</span>
<span class="hl str">			ibd0_cM_m2 = (hypo_cM_length_map2 - (cm_length_m2_IBD1 + cm_length_m2_IBD2)) / hypo_cM_length_map2,</span>
<span class="hl str">			ibd1_cM_m2 = cm_length_m2_IBD1 / hypo_cM_length_map2,</span>
<span class="hl str">			ibd2_cM_m2 = cm_length_m2_IBD2 / hypo_cM_length_map2)</span>
<span class="hl str"></span>
<span class="hl str">	cM_treshold &lt;- 0.2</span>
<span class="hl str">	summary_filterd &lt;- converted_segments %&gt;%</span>
<span class="hl str">	filter(length_cM_m1 &gt; cM_treshold &amp; length_cM_m2 &gt; cM_treshold) %&gt;% </span>
<span class="hl str">	ungroup() %&gt;%</span>
<span class="hl str">	group_by(PAIR, TYPE) %&gt;%</span>
<span class="hl str">	summarise(seq_length = sum(length_bp),</span>
<span class="hl str">				n_mark = sum(NMARKERS),</span>
<span class="hl str">				cm_length_m1 = sum(length_cM_m1),</span>
<span class="hl str">				cm_length_m2 = sum(length_cM_m2)) %&gt;%</span>
<span class="hl str">	ungroup() %&gt;%</span>
<span class="hl str">	pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:cm_length_m2, values_fill = 0) %&gt;% </span>
<span class="hl str">	left_join(segments_summary, .,  ) %&gt;% </span>
<span class="hl str">	mutate(across(.cols = seq_length_IBD1:cm_length_m2_IBD2, .fns = na_to_zero)) %&gt;% </span>
<span class="hl str">	mutate(IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,</span>
<span class="hl str">			IBD1_manual = n_mark_IBD1 / NMARK,</span>
<span class="hl str">			IBD2_manual = n_mark_IBD2 / NMARK,</span>
<span class="hl str">			icheck_0 = IBD0_manual - IBD0,</span>
<span class="hl str">			icheck_1 = IBD1_manual - IBD1,</span>
<span class="hl str">			icheck_2 = IBD2 - IBD2_manual,</span>
<span class="hl str">			# compile ibd by sequence map</span>
<span class="hl str">			ibd0_bp = (hypo_bp_length - (seq_length_IBD1 + seq_length_IBD2)) / hypo_bp_length,</span>
<span class="hl str">			ibd1_bp = seq_length_IBD1 / hypo_bp_length,</span>
<span class="hl str">			ibd2_bp = seq_length_IBD2 / hypo_bp_length,</span>
<span class="hl str">			# compile ibd by genetic map 1</span>
<span class="hl str">			ibd0_cM_m1 = (hypo_cM_length_map1 - (cm_length_m1_IBD1 + cm_length_m1_IBD2)) / hypo_cM_length_map1,</span>
<span class="hl str">			ibd1_cM_m1 = cm_length_m1_IBD1 / hypo_cM_length_map1,</span>
<span class="hl str">			ibd2_cM_m1 = cm_length_m1_IBD2 / hypo_cM_length_map1,</span>
<span class="hl str">			# compile ibd by genetic map 2</span>
<span class="hl str">			ibd0_cM_m2 = (hypo_cM_length_map2 - (cm_length_m2_IBD1 + cm_length_m2_IBD2)) / hypo_cM_length_map2,</span>
<span class="hl str">			ibd1_cM_m2 = cm_length_m2_IBD1 / hypo_cM_length_map2,</span>
<span class="hl str">			ibd2_cM_m2 = cm_length_m2_IBD2 / hypo_cM_length_map2)</span>
<span class="hl str"></span>
<span class="hl str">	write_tsv(x = converted_segments, file = truffle_conv)</span>
<span class="hl str">	write_tsv(x = control, file = truffle_sum)</span>
<span class="hl str">	write_tsv(x = summary_filterd, file = truffle_filt)</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
