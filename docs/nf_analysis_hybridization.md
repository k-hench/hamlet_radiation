---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 9) Analysis VII (hybridization)

## Summary

The population recombination rate is estimated within the [**nextflow**](https://www.nextflow.io/) script `analysis_hybridization.nf` (located under `$BASE_DIR/nf/09_analysis_hybridization/`), which runs on the XX data set.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:500px; margin:auto;">

</div>

## Details of `analysis_hybridization.nf`

### Data preparation

The nextflow script starts by opening the genotype data.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 9.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.into<span class="hl opt">{</span> vcf_loc1<span class="hl opt">;</span> vcf_loc2<span class="hl opt">;</span> vcf_loc3 <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.2</span>
<span class="hl slc">// initialize location channel</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">(</span> <span class="hl str">&quot;bel&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;hon&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;pan&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> locations_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.3</span>
<span class="hl slc">// define location specific sepcies set</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;ind&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;may&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> bel_spec1_ch<span class="hl opt">;</span> bel_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;abe&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;gum&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">4</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">5</span><span class="hl opt">,</span> <span class="hl str">&quot;ran&quot;</span><span class="hl opt">], [</span><span class="hl num">6</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> hon_spec1_ch<span class="hl opt">;</span> hon_spec2_ch <span class="hl opt">}</span>
<span class="hl kwa">Channel</span>.from<span class="hl opt">( [[</span><span class="hl num">1</span><span class="hl opt">,</span> <span class="hl str">&quot;nig&quot;</span><span class="hl opt">], [</span><span class="hl num">2</span><span class="hl opt">,</span> <span class="hl str">&quot;pue&quot;</span><span class="hl opt">], [</span><span class="hl num">3</span><span class="hl opt">,</span> <span class="hl str">&quot;uni&quot;</span><span class="hl opt">]] )</span>.into<span class="hl opt">{</span> pan_spec1_ch<span class="hl opt">;</span> pan_spec2_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.4</span>
<span class="hl slc">// prepare pairwise new_hybrids</span>
<span class="hl slc">// ------------------------------</span>
<span class="hl com">/* (create all possible species pairs depending on location</span>
<span class="hl com">   and combine with genotype subset (for the respective location))*/</span>
bel_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;bel&quot;</span> <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span> vcf_loc1 <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>bel_spec1_ch<span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>bel_spec2_ch<span class="hl opt">)</span>
    .filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">3</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] }</span>
    .map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">]}</span>
hon_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;hon&quot;</span> <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span> vcf_loc2 <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>hon_spec1_ch<span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>hon_spec2_ch<span class="hl opt">)</span>
    .filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">3</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] }</span>
    .map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">]}</span>
pan_pairs_ch <span class="hl opt">=</span> <span class="hl kwa">Channel</span>.from<span class="hl opt">(</span> <span class="hl str">&quot;pan&quot;</span> <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span> vcf_loc3 <span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>pan_spec1_ch<span class="hl opt">)</span>
    .combine<span class="hl opt">(</span>pan_spec2_ch<span class="hl opt">)</span>
    .filter<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">3</span><span class="hl opt">] &lt;</span> it<span class="hl opt">[</span><span class="hl num">5</span><span class="hl opt">] }</span>
    .map<span class="hl opt">{</span> it<span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">,</span><span class="hl num">1</span><span class="hl opt">,</span><span class="hl num">2</span><span class="hl opt">,</span><span class="hl num">4</span><span class="hl opt">,</span><span class="hl num">6</span><span class="hl opt">]}</span>
bel_pairs_ch.concat<span class="hl opt">(</span> hon_pairs_ch<span class="hl opt">,</span> pan_pairs_ch  <span class="hl opt">)</span>.set <span class="hl opt">{</span> all_fst_pairs_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.5</span>
<span class="hl slc">// comute pairwise fsts for SNP filtering</span>
<span class="hl kwa">process</span> fst_run <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g45m_fst_run&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">)</span> <span class="hl kwa">from</span> all_fst_pairs_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${vcf[0]}</span><span class="hl str">&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.fst.tsv.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">.pop&quot;</span><span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.pop&quot;</span><span class="hl opt">)</span> <span class="hl kwa">into</span> fst_SNPS

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| grep</span> <span class="hl ipl">${spec1}${loc}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">.pop</span>
<span class="hl str">	vcfsamplenames</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">| grep</span> <span class="hl ipl">${spec2}${loc}</span> <span class="hl str">&gt;</span> <span class="hl ipl">${spec2}${loc}</span><span class="hl str">.pop</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		 --weir-fst-pop</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">.pop \</span>
<span class="hl str">		 --weir-fst-pop</span> <span class="hl ipl">${spec2}${loc}</span><span class="hl str">.pop \</span>
<span class="hl str">		 --stdout | gzip &gt;</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.fst.tsv.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.6</span>
<span class="hl slc">// select the 800 most differentiated SNPs for each population pair</span>
<span class="hl kwa">process</span> filter_fst <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_8g15m_filter_fst&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> fst <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop1 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop2 <span class="hl opt">)</span> <span class="hl kwa">from</span> fst_SNPS

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop1 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*SNPs.snps&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> filter_SNPs


	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/filter_snps.R</span> <span class="hl ipl">${fst}</span> <span class="hl str">800</span> <span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.7</span>
<span class="hl slc">// filter the SNP set by min distance (5kb), than randomly pick 80 SNPs</span>
<span class="hl slc">// then reformat newhybrid input</span>
<span class="hl kwa">process</span> prep_nh_input <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_8g15m_prep_nh&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop1 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pop2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> snps <span class="hl opt">)</span> <span class="hl kwa">from</span> filter_SNPs


	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_individuals.txt&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.80SNPs.txt&quot;</span><span class="hl opt">)</span>  <span class="hl kwa">into</span> newhybrids_input


	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">  --gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">	--keep</span> <span class="hl ipl">${pop1}</span> <span class="hl str">\</span>
<span class="hl str">	--keep</span> <span class="hl ipl">${pop2}</span> <span class="hl str">\</span>
<span class="hl str">	--thin 5000 \</span>
<span class="hl str">	--out newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span> <span class="hl str">\</span>
<span class="hl str">	--positions</span> <span class="hl ipl">${snps}</span> <span class="hl str">\</span>
<span class="hl str">	--recode</span>
<span class="hl str"></span>
<span class="hl str">	grep &#39;#&#39; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.recode.vcf &gt; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs.vcf</span>
<span class="hl str">	grep -v &#39;#&#39; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.recode.vcf | \</span>
<span class="hl str">		shuf -n 80 | \</span>
<span class="hl str">		sort -k 1 -k2 &gt;&gt; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs.vcf</span>
<span class="hl str"></span>
<span class="hl str">	grep &#39;#CHROM&#39; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs.vcf | \</span>
<span class="hl str">		cut -f 10- | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">t/</span><span class="hl esc">\\</span><span class="hl str">n/g&#39; &gt; newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs_individuals.txt</span>
<span class="hl str"></span>
<span class="hl str">	/usr/bin/java -Xmx1024m -Xms512M \</span>
<span class="hl str">		-jar \$SFTWR/PGDSpider/PGDSpider2-cli.jar \</span>
<span class="hl str">		-inputfile newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs.vcf \</span>
<span class="hl str">		-inputformat VCF \</span>
<span class="hl str">		-outputfile newHyb.</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">.80SNPs.txt \</span>
<span class="hl str">		-outputformat NEWHYBRIDS \</span>
<span class="hl str">		-spid \$BASE_DIR/ressources/vcf2nh.spid</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 9.8</span>
<span class="hl slc">// Run new hybrids</span>
<span class="hl slc">// (copy of nh_input is needed because nh can&#39;t read links)</span>
<span class="hl kwa">process</span> run_nh <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g15h4x_run_nh&#39;</span>
	<span class="hl kwb">tag</span> <span class="hl str">&quot;</span><span class="hl ipl">${spec1}${loc}</span><span class="hl str">-</span><span class="hl ipl">${spec2}${loc}</span><span class="hl str">&quot;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/newhyb/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> loc <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec1 <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> spec2 <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> inds <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> snps <span class="hl opt">)</span> <span class="hl kwa">from</span> newhybrids_input

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;nh_input/NH.Results/newHyb.*/*_individuals.txt&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;nh_input/NH.Results/newHyb.*/*_PofZ.txt&quot;</span> <span class="hl opt">)</span>  <span class="hl kwa">into</span> newhybrids_output

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	mkdir -p nh_input</span>
<span class="hl str">	cp</span> <span class="hl ipl">${snps}</span> <span class="hl str">nh_input/</span><span class="hl ipl">${snps}</span>
<span class="hl str">	cp</span> <span class="hl ipl">${inds}</span> <span class="hl str">nh_input/</span><span class="hl ipl">${inds}</span>
<span class="hl str"></span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/run_newhybrids.R</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
</div>

Finally, we are done with XX.

---
