---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 11) Analysis IX (Allele Age)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/11_analysis_allele_age
source ../sh/nextflow_alias.sh
nf_run_aa
```

## Summary

The allele age is estimated within the [**nextflow**](https://www.nextflow.io/) script `analysis_allele_age.nf` (located under `$BASE_DIR/nf/11_analysis_allele_age/`) which runs on the _SNPs only_ data set.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:800px; margin:auto;">

</div>

## Details of `analysis_allele_age.nf`

### Setup

The nextflow script starts by opening the genotype data.

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 11.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

As the data is going to be split by linkage group, we create a channel for the individual LGs.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.2</span>
<span class="hl slc">// initialize LGs</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;01&#39;</span>..<span class="hl str">&#39;09&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;10&#39;</span>..<span class="hl str">&#39;19&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;20&#39;</span>..<span class="hl str">&#39;24&#39;</span><span class="hl opt">) )</span>
	.map<span class="hl opt">{</span><span class="hl str">&quot;LG&quot;</span> <span class="hl opt">+</span> it<span class="hl opt">}</span>
	.combine<span class="hl opt">(</span> vcf_ch <span class="hl opt">)</span>
	.set<span class="hl opt">{</span> lg_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Several things are happening in the next step:

- First the genotypes are split by LG and a new info field is added to the `vcf` file to store information about the ancestral state of each SNP.
- The second step checks for each SNP if it is invariant across all Serranid outgroup samples (in that case the outgroup allele is considered ancestral).
- Then, allele frequencies are computed as a fallback clue - if a SNP is variant in the outgroup, the major allele is then set as ancestral allele.
- Lastly, the ancestral state information is added to the genotype file.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.3</span>
<span class="hl slc">// subset the genotypes by LG</span>
<span class="hl slc">// and add ancestral allele annotation</span>
<span class="hl kwa">process</span> prepare_vcf <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_20g2h_prepare_vcf&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> vcfidx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> lg_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${lg}</span><span class="hl str">_integer.vcf&quot;</span> <span class="hl opt">)</span>  <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_prep_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# subset by LG and add AC info field</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">	  --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		--chr</span> <span class="hl ipl">${lg}</span> <span class="hl str">\</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">(##FORMAT=&lt;ID=GT,Number=1,Type=String,Description=&quot;Phased Genotype&quot;&gt;</span><span class="hl esc">\\</span><span class="hl str">)/</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">n##INFO=&lt;ID=AC,Number=A,Type=Integer,Description=&quot;Allele count in genotypes, for each ALT allele, in the same order as listed&quot;&gt;/&#39; | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	# determine ancestral state based on invariant sites in outgoup</span>
<span class="hl str">	zcat</span> <span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz | \</span>
<span class="hl str">		grep -v &quot;^#&quot; | \</span>
<span class="hl str">		awk -v OFS=&quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; \</span>
<span class="hl str">		&#39;BEGIN{print &quot;#CHROM&quot;,&quot;FROM&quot;,&quot;TO&quot;,&quot;AA&quot;}</span>
<span class="hl str">		{o1=substr(\$107, 1, 1);</span>
<span class="hl str">		o2=substr(\$107, 3, 3);</span>
<span class="hl str">		o3=substr(\$167, 1, 1);</span>
<span class="hl str">		o4=substr(\$167, 3, 3);</span>
<span class="hl str">		o5=substr(\$179, 1, 1);</span>
<span class="hl str">		o6=substr(\$179, 3, 3);</span>
<span class="hl str">		if (o1 == o2 &amp;&amp; o3 == o4 &amp;&amp; o5 == o6 &amp;&amp; o1 == o3 &amp;&amp; o1 == o5){</span>
<span class="hl str">		aa = \$(4+o1)} else {aa = &quot;.&quot;};</span>
<span class="hl str">		print \$1,\$2,\$2,aa}&#39; &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_annotations.bed</span>
<span class="hl str"></span>
<span class="hl str">	# determine allele frquencies</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz \</span>
<span class="hl str">		--freq \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		sed &#39;s/{ALLELE:FREQ}/ALLELE1</span><span class="hl esc">\\</span><span class="hl str">tALLELE2/&#39; &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_allele_counts.tsv</span>
<span class="hl str"></span>
<span class="hl str">	# determine ancestral state for variant sites in outgoup based on allele freq</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/major_allele.R</span> <span class="hl ipl">${lg}</span><span class="hl str">_allele_counts.tsv</span> <span class="hl ipl">${lg}</span><span class="hl str">_annotations.bed</span>
<span class="hl str"></span>
<span class="hl str">	bgzip</span> <span class="hl ipl">${lg}</span><span class="hl str">_annotations_maj.bed</span>
<span class="hl str"></span>
<span class="hl str">	tabix -s 1 -b 2 -e 3</span> <span class="hl ipl">${lg}</span><span class="hl str">_annotations_maj.bed.gz</span>
<span class="hl str"></span>
<span class="hl str">	# add ancestral state annotation</span>
<span class="hl str">	zcat</span> <span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz | \</span>
<span class="hl str">		vcf-annotate -a</span> <span class="hl ipl">${lg}</span><span class="hl str">_annotations_maj.bed.gz \</span>
<span class="hl str">		-d key=INFO,ID=AA,Number=1,Type=String,Description=&#39;Ancestral Allele&#39; \</span>
<span class="hl str">		-c CHROM,FROM,TO,INFO/AA | \</span>
<span class="hl str">		sed &#39;s/LG//g&#39;  \</span>
<span class="hl str">		&gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_integer.vcf</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Based on the ancestral state information, the genotypes are recoded such that the ancestral allele is set as the new reference allele of the `vcf` file.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.4</span>
<span class="hl slc">// re-write ancestral state in vcf</span>
<span class="hl kwa">process</span> set_ancestral_states <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_2g15m_ancestral_states&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../1_genotyping/5_ancestral_allele&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> <span class="hl opt">(</span> vcf_prep_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${lg}</span><span class="hl str">_aa.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_aa_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	java -jar \$SFTWR/jvarkit/dist/vcffilterjdk.jar \</span>
<span class="hl str">		-f \$BASE_DIR/js/script.js</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_aa.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

After this, the outgroup samples are removed from the data and the remaining data set is filtered to remove sites that are invariant within the hamlets.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.5</span>
<span class="hl slc">// filter vcf to remove invariant sites in hamlets</span>
<span class="hl kwa">process</span> create_positions <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_20g2h_create_positions&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/sliding_phylo/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> <span class="hl opt">(</span> vcf_aa_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${lg}</span><span class="hl str">_aa_h_variant.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;</span><span class="hl ipl">${lg}</span><span class="hl str">_positions.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> positions_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo -e &quot;20478tabhon</span><span class="hl esc">\\</span><span class="hl str">n28393torpan</span><span class="hl esc">\\</span><span class="hl str">ns_tort_3torpan&quot; &gt; outgr.pop</span>
<span class="hl str"></span>
<span class="hl str">	# keeping only sites that are variant within hamlets</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--remove outgr.pop \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		vcftools \</span>
<span class="hl str">		--gzvcf - \</span>
<span class="hl str">		--mac 1 \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_aa_no_outgroup.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	zcat</span> <span class="hl ipl">${lg}</span><span class="hl str">_aa_no_outgroup.vcf.gz | \</span>
<span class="hl str">		grep -v &quot;^#&quot; | \</span>
<span class="hl str">		cut -f 1,2 | \</span>
<span class="hl str">		head -n -1 &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_positions_prep.txt</span>
<span class="hl str"></span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--positions</span> <span class="hl ipl">${lg}</span><span class="hl str">_positions_prep.txt \</span>
<span class="hl str">		--recode \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		bgzip &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_aa_h_variant.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	cut -f 2</span> <span class="hl ipl">${lg}</span><span class="hl str">_positions_prep.txt &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_positions.txt</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Since a single `GEVA` run can take quite some time, the data is split further into chunks of 25k SNPs 


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.5</span>
<span class="hl slc">// prepare the age estimation by splitting the vcf</span>
<span class="hl slc">// (all in one takes too long...)</span>
<span class="hl kwa">process</span> pre_split <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_2g2h_pre_split&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/geva/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pos <span class="hl opt">)</span> <span class="hl kwa">from</span> <span class="hl opt">(</span> positions_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;pre_positions/pre_*&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.bin&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.marker.txt&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.sample.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> geva_setup_ch <span class="hl opt">)</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;inner_pos.txt&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> ccf_vcf_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	mkdir -p pre_positions</span>
<span class="hl str"></span>
<span class="hl str">	head -n -1</span> <span class="hl ipl">${pos}</span> <span class="hl str">| \</span>
<span class="hl str">	 tail -n +2  &gt; inner_pos.txt</span>
<span class="hl str"></span>
<span class="hl str">	split inner_pos.txt -a 4 -l 25000 -d pre_positions/pre_</span>
<span class="hl str"></span>
<span class="hl str">	r=\$(awk -v k=</span><span class="hl ipl">${lg}</span> <span class="hl str">&#39;\$1 == k {print \$4}&#39; \$BASE_DIR/ressources/avg_rho_by_LG.tsv)</span>
<span class="hl str"></span>
<span class="hl str">	geva_v1beta \</span>
<span class="hl str">		--vcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">--rec \$r --out</span> <span class="hl ipl">${lg}</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then `GEVA` is run on those genotype subsets to actually estimate the allele ages.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.6</span>
<span class="hl slc">// run geva on vcf subsets</span>
<span class="hl kwa">process</span> run_geva <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_30g15h6x_run_geva&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pos <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> bin <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> marker <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> sample <span class="hl opt">)</span> <span class="hl kwa">from</span> geva_setup_ch.transpose<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.sites.txt.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.pairs.txt.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> output_split_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">  pref=\$(echo &quot;</span><span class="hl ipl">${pos}</span><span class="hl str">&quot; | sed &#39;s=^.*/A==; s=pre_positions/pre_==&#39;)</span>
<span class="hl str"></span>
<span class="hl str">	mkdir -p sub_positions sub_results</span>
<span class="hl str"></span>
<span class="hl str">	split</span> <span class="hl ipl">${pos}</span> <span class="hl str">-a 4 -l 250 -d sub_positions/sub_pos_\</span><span class="hl ipl">${pref}</span><span class="hl str">_</span>
<span class="hl str"></span>
<span class="hl str">	r=\$(awk -v k=</span><span class="hl ipl">${lg}</span> <span class="hl str">&#39;\$1 == k {print \$4}&#39; \$BASE_DIR/ressources/avg_rho_by_LG.tsv)</span>
<span class="hl str"></span>
<span class="hl str">	for sp in \$(ls sub_positions/sub_pos_\</span><span class="hl ipl">${pref}</span><span class="hl str">_*); do</span>
<span class="hl str">		run_id=\$(echo \$sp | sed &quot;s=sub_positions/sub_pos_\</span><span class="hl ipl">${pref}</span><span class="hl str">_==&quot;)</span>
<span class="hl str"></span>
<span class="hl str">		geva_v1beta \</span>
<span class="hl str">			 -t 6 \</span>
<span class="hl str">			 -i</span> <span class="hl ipl">${bin}</span> <span class="hl str">\</span>
<span class="hl str">			 -o sub_results/</span><span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">_\</span><span class="hl ipl">${run_id}</span><span class="hl str">\</span>
<span class="hl str">			 --positions \$sp \</span>
<span class="hl str">			 --Ne 30000 \</span>
<span class="hl str">			 --mut 3.7e-08 \</span>
<span class="hl str">			 --hmm \$SFTWR/geva/hmm/hmm_initial_probs.txt \$SFTWR/geva/hmm/hmm_emission_probs.txt</span>
<span class="hl str"></span>
<span class="hl str">		tail -n +2 sub_results/</span><span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">_\</span><span class="hl ipl">${run_id}</span><span class="hl str">.sites.txt &gt;&gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">.sites.txt</span>
<span class="hl str">		tail -n +2 sub_results/</span><span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">_\</span><span class="hl ipl">${run_id}</span><span class="hl str">.pairs.txt &gt;&gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">.pairs.txt</span>
<span class="hl str">	done</span>
<span class="hl str"></span>
<span class="hl str">	gzip</span> <span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">.sites.txt</span>
<span class="hl str">	gzip</span> <span class="hl ipl">${lg}</span><span class="hl str">_\</span><span class="hl ipl">${pref}</span><span class="hl str">.pairs.txt</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

And finally, the results of the separate chunks are gathered and compiled into a single output file.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 11.7</span>
<span class="hl slc">// collect results by lg</span>
<span class="hl kwa">process</span> collect_by_lg <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_2g2h_collect&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/geva/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> sites <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> pairs <span class="hl opt">)</span> <span class="hl kwa">from</span> output_split_ch.groupTuple<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.sites.txt.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.pairs.txt.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> output_lg_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	echo &quot;MarkerID Clock Filtered N_Concordant N_Discordant PostMean PostMode PostMedian&quot; &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">.sites.txt</span>
<span class="hl str">	echo &quot;MarkerID Clock SampleID0 Chr0 SampleID1 Chr1 Shared Pass SegmentLHS SegmentRHS Shape Rate&quot; &gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">.pairs.txt</span>
<span class="hl str"></span>
<span class="hl str">	zcat</span> <span class="hl ipl">${sites}</span>  <span class="hl str">&gt;&gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">.sites.txt</span>
<span class="hl str">	zcat</span> <span class="hl ipl">${pairs}</span>  <span class="hl str">&gt;&gt;</span> <span class="hl ipl">${lg}</span><span class="hl str">.pairs.txt</span>
<span class="hl str"></span>
<span class="hl str">	gzip</span> <span class="hl ipl">${lg}</span><span class="hl str">.sites.txt</span>
<span class="hl str">	gzip</span> <span class="hl ipl">${lg}</span><span class="hl str">.pairs.txt</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
