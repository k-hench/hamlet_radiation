---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 15) Analysis XIII (diversity with and without outlier regions)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/15_analysis_pi
source ../../sh/nextflow_alias.sh
nf_run_pi
```

## Summary

The diversity is computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_pi.nf` (located under `$BASE_DIR/nf/15_analysis_pi/`).
It takes the _all BP_ data set and computes $\pi$.

## Details of `analysis_pi.nf`

### Setup

The nextflow script starts by opening the genotype data and feeding it into a stream.

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// This pipeline includes the analysis run on the</span>
<span class="hl slc">//   all callable sites data sheet (pi).</span>

<span class="hl slc">// git 15.1</span>
<span class="hl slc">// load genotypes</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/filterd.allBP.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_pi_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

To split the analysis by linkage group, we initialize a LG-channel.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 15.2</span>
<span class="hl slc">// initialize LGs</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;01&#39;</span>..<span class="hl str">&#39;09&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;10&#39;</span>..<span class="hl str">&#39;19&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;20&#39;</span>..<span class="hl str">&#39;24&#39;</span><span class="hl opt">) )</span>
	.set<span class="hl opt">{</span> lg_pi_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

We also want to compute the diversity twice, once *as is* and once without the outlier regions.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 15.3</span>
<span class="hl slc">// set outlier window mode</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;_no_outlier&quot;</span> <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> outlier_mode_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Furthermore, we want to analyze diversity at two scales (10kb and 50kb).


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 15.4</span>
<span class="hl slc">// init slining window resolutions</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl num">1</span><span class="hl opt">,</span> <span class="hl num">5</span> <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> kb_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Then we filter the genotypes, depending on population, linkage group and outlier mode.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 15.5</span>
<span class="hl slc">// init all sampled populations (for pi)</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span><span class="hl str">&#39;indbel&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;maybel&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;nigbel&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;puebel&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;unibel&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;abehon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;gumhon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;nighon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;puehon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;ranhon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;unihon&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;nigpan&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;puepan&#39;</span><span class="hl opt">,</span> <span class="hl str">&#39;unipan&#39;</span><span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> vcf_pi_ch <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> kb_ch <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> lg_pi_ch <span class="hl opt">)</span>
	.combine<span class="hl opt">(</span> outlier_mode_ch <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> input_ch <span class="hl opt">}</span>

<span class="hl slc">// filter genotypes and convert formats</span>
<span class="hl kwa">process</span> recode_genotypes <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g15h_recode&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlier <span class="hl opt">)</span> <span class="hl kwa">from</span> input_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlier <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.geno.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;pop.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> geno_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	if [ &quot;</span><span class="hl ipl">${outlier}</span><span class="hl str">&quot; == &quot;_no_outlier&quot; ];then</span>
<span class="hl str">		tail -n +2 \$BASE_DIR/2_analysis/summaries/fst_outliers_998.tsv | \</span>
<span class="hl str">			cut -f 2,3,4 &gt; outlier.bed </span>
<span class="hl str">		SUBSET=&quot;--exclude-bed outlier.bed&quot;</span>
<span class="hl str">	else</span>
<span class="hl str">		SUBSET=&quot;&quot;</span>
<span class="hl str">	fi</span>
<span class="hl str"></span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| \</span>
<span class="hl str">		grep</span> <span class="hl ipl">${spec}</span> <span class="hl str">&gt; pop.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--keep pop.txt \</span>
<span class="hl str">		--chr LG</span><span class="hl ipl">${lg}</span> <span class="hl str">\</span>
<span class="hl str">		\$SUBSET \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | bgzip &gt;</span> <span class="hl ipl">${spec}${outlier}</span><span class="hl str">.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \</span>
<span class="hl str">		-i</span> <span class="hl ipl">${spec}${outlier}</span><span class="hl str">.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz | \</span>
<span class="hl str">		gzip &gt;</span> <span class="hl ipl">${spec}${outlier}</span><span class="hl str">.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.geno.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then we calculate pi for each population.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 15.6</span>
<span class="hl slc">// calculate pi per species</span>
<span class="hl kwa">process</span> pi_per_spec <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g15h_pi&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlier <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> geno <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop <span class="hl opt">)</span> <span class="hl kwa">from</span> geno_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec}${outlier}</span><span class="hl str">.</span><span class="hl ipl">${kb}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*0kb.csv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> pi_lg_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	awk &#39;{print \$1&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot;substr(\$1,length(\$1)-5,length(\$1)-1)}&#39;</span> <span class="hl ipl">${pop}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${spec}</span><span class="hl str">.pop</span>
<span class="hl str"></span>
<span class="hl str">	python \$SFTWR/genomics_general/popgenWindows.py \</span>
<span class="hl str">		-w</span> <span class="hl ipl">${kb}</span><span class="hl str">0000 \</span>
<span class="hl str">		-s</span> <span class="hl ipl">${kb}</span><span class="hl str">000 \</span>
<span class="hl str">		--popsFile</span> <span class="hl ipl">${spec}</span><span class="hl str">.pop \</span>
<span class="hl str">		-p</span> <span class="hl ipl">${spec}</span> <span class="hl str">\</span>
<span class="hl str">		-g</span> <span class="hl ipl">${geno}</span> <span class="hl str">\</span>
<span class="hl str">		-o pi.</span><span class="hl ipl">${spec}${outlier}</span><span class="hl str">.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.</span><span class="hl ipl">${kb}</span><span class="hl str">0kb.csv.gz \</span>
<span class="hl str">		-f phased \</span>
<span class="hl str">		--writeFailedWindows \</span>
<span class="hl str">		-T 1</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Finally, we merge the output from the individual linkage groups.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl kwa">process</span> merge_pi <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_32g4h_pi_merge&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec_outlier_kb}</span><span class="hl str">&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/pi/</span><span class="hl ipl">${kb[0]}</span><span class="hl str">0k&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec_outlier_kb <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> kb <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pi <span class="hl opt">)</span> <span class="hl kwa">from</span> pi_lg_ch.groupTuple<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;pi.</span><span class="hl ipl">${spec_outlier_kb}</span><span class="hl str">0kb.tsv.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> pi_output_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo -e &quot;CHROM</span><span class="hl esc">\\</span><span class="hl str">tBIN_START</span><span class="hl esc">\\</span><span class="hl str">tBIN_END</span><span class="hl esc">\\</span><span class="hl str">tBIN_MID_SITES</span><span class="hl esc">\\</span><span class="hl str">tN_SITES</span><span class="hl esc">\\</span><span class="hl str">tPI&quot; &gt; pi.</span><span class="hl ipl">${spec_outlier_kb}</span><span class="hl str">0kb.tsv</span>
<span class="hl str"></span>
<span class="hl str">	for k in \$(ls *.csv.gz); do </span>
<span class="hl str">		zcat \$k | \</span>
<span class="hl str">			grep -v &quot;scaff&quot; | \</span>
<span class="hl str">			sed -s &quot;s/,/</span><span class="hl esc">\t</span><span class="hl str">/g&quot;  &gt;&gt; pi.</span><span class="hl ipl">${spec_outlier_kb}</span><span class="hl str">0kb.tsv</span>
<span class="hl str">	done</span>
<span class="hl str">	</span>
<span class="hl str">	gzip pi.</span><span class="hl ipl">${spec_outlier_kb}</span><span class="hl str">0kb.tsv</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
