---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---



# (git 13) Analysis  XI (Whole Genome Phylogenies)

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

<span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// this pipeline was not actually run using nexflow,</span>
<span class="hl slc">// but managed manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>

<span class="hl slc">// Hamlet phylogeny</span>
<span class="hl slc">// ----------------</span>

<span class="hl slc">// git 13.1</span>
<span class="hl slc">// open the SNP data set</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.into<span class="hl opt">{</span> vcf_hypo_whg_ch<span class="hl opt">;</span> vcf_serr_whg_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Next, also a file containing the sample IDs excluding the samples identified as hybrids in git 9 is loaded.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// RAxML analysis, Serranus-rooted</span>
<span class="hl slc">// -------------------------------</span>
<span class="hl slc">// git 13.2</span>
<span class="hl slc">// open the sample-list (excluding hybrid samples)</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/samples_hybrids.txt&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> hybrids_file <span class="hl opt">}</span>
</code>
</pre>
</div>

As a preparation for running `raxml`, the genotype file is subset to exclude the hybrids.
Then, heterozygous sites are masked and the data is filtered for minimal allele count and physical distance thresholds.
Finally, the genotypes a indirectly converted to `fasta` format.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.3</span>
<span class="hl slc">// subset data and convert to fasta for raxml</span>
<span class="hl kwa">process</span> serr_whg_genotypes <span class="hl opt">{</span>
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> hybrids <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_serr_whg_ch.combine<span class="hl opt">(</span> hybrids_file <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyS_n_0.33_mac4_5kb.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> raxml_serr_genotypes_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# Remove hybrids from genotype data (SNPs only)</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">	  --gzvcf</span>  <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">	  --remove</span> <span class="hl ipl">${hybrids}</span> <span class="hl str">\</span>
<span class="hl str">	  --recode \</span>
<span class="hl str">	  --stdout | \</span>
<span class="hl str">	  gzip &gt; hyS.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	# Mask heterozygous genotypes as unknown</span>
<span class="hl str">	zcat &lt; hyS.vcf.gz | \</span>
<span class="hl str">	  sed -e s/&quot;1|0&quot;/&quot;.|.&quot;/g -e s/&quot;0|1&quot;/&quot;.|.&quot;/g | \</span>
<span class="hl str">	  gzip &gt; hyS_n.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	# Apply missingness, allele count and distance filters</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">	  --gzvcf hyS_n.vcf.gz \</span>
<span class="hl str">	  --max-missing 0.33 \</span>
<span class="hl str">	  --mac 4 \</span>
<span class="hl str">	  --thin 5000 \</span>
<span class="hl str">	  --recode \</span>
<span class="hl str">	  --out hyS_n_0.33_mac4_5kb</span>
<span class="hl str"></span>
<span class="hl str">	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)</span>
<span class="hl str">	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl</span>
<span class="hl str"></span>
<span class="hl str">	vcf-to-tab &lt; hyS_n_0.33_mac4_5kb.vcf &gt; hyS_n_0.33_mac4_5kb.tab</span>
<span class="hl str">	</span>
<span class="hl str">	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i hyS_n_0.33_mac4_5kb.tab &gt; hyS_n_0.33_mac4_5kb.fas</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then, `raxml` can be run on the `fasta` formated genotypes.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.4</span>
<span class="hl slc">// run raxml (Serranus-rooted)</span>
<span class="hl kwa">process</span> serr_whg_raxml <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/raxml/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> fas <span class="hl opt">)</span> <span class="hl kwa">from</span> raxml_serr_genotypes_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyS_n_0.33_mac4_5kb.raxml.support&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> raxml_serr_whg_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# Reconstruct phylogeny</span>
<span class="hl str">	# Note: number of invariant sites for Felsenstein correction was calculated as number of</span>
<span class="hl str">	# variant sites in alignment (109,660) / genome-wide proportion of variant sites</span>
<span class="hl str">	# (0.05) * genome-wide proportion of invariant sites (0.95)</span>
<span class="hl str">	raxml-NG --all \</span>
<span class="hl str">	  --msa hyS_n_0.33_mac4_5kb.fas \</span>
<span class="hl str">	  --model GTR+G+ASC_FELS{2083540} \</span>
<span class="hl str">	  --tree pars{20},rand{20} \</span>
<span class="hl str">	  --bs-trees 100 \</span>
<span class="hl str">	  --threads 24 \</span>
<span class="hl str">	  --worker 4 \</span>
<span class="hl str">	  --seed 123 \</span>
<span class="hl str">	  --prefix hyS_n_0.33_mac4_5kb</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

---

The same general aproach is used for the phylogeny excluding the Serranus outgroup samples.
For this, a different sample list (also excluding the outgroup samples) is loaded.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// RAxML analysis, floridae-rooted</span>
<span class="hl slc">// -------------------------------</span>
<span class="hl slc">// git 13.5</span>
<span class="hl slc">// open the sample-list (excluding hybrid and Serranus samples)</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/samples_155.txt&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> hamlet_file <span class="hl opt">}</span>
</code>
</pre>
</div>

Like in git 13.3, the genotypes are subset and converted to `fasta` format.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.6</span>
<span class="hl slc">// subset data and convert to fasta for raxml</span>
<span class="hl kwa">process</span> hypo_whg_genotypes <span class="hl opt">{</span>
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> hamlets <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_hypo_whg_ch.combine<span class="hl opt">(</span>hamlet_file<span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp155_n_0.33_mac4_5kb.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> raxml_hypo_genotypes_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# Remove hybrid and Serranus samples from genotype data (SNPs only)</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">	  --gzvcf</span> <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">	  --remove</span> <span class="hl ipl">${hamlets}</span> <span class="hl str">\</span>
<span class="hl str">	  --recode \</span>
<span class="hl str">	  --stdout | \</span>
<span class="hl str">	  gzip &gt; hyp155.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	# Mask heterozygous genotypes as unknown</span>
<span class="hl str">	zcat &lt; hyp155.vcf.gz | \</span>
<span class="hl str">	  sed -e s/&quot;1|0&quot;/&quot;.|.&quot;/g -e s/&quot;0|1&quot;/&quot;.|.&quot;/g | \</span>
<span class="hl str">	  gzip &gt; hyp155_n.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	# Apply missingness, allele count and distance filters</span>
<span class="hl str">	vcftools \</span>
<span class="hl str">	  --gzvcf hyp155_n.vcf.gz \</span>
<span class="hl str">	  --max-missing 0.33 \</span>
<span class="hl str">	  --mac 4 \</span>
<span class="hl str">	  --thin 5000 \</span>
<span class="hl str">	  --recode \</span>
<span class="hl str">	  --out hyp155_n_0.33_mac4_5kb</span>
<span class="hl str"></span>
<span class="hl str">	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)</span>
<span class="hl str">	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl</span>
<span class="hl str"></span>
<span class="hl str">	vcf-to-tab &lt; hyp155_n_0.33_mac4_5kb.vcf &gt; hyp155_n_0.33_mac4_5kb.tab</span>
<span class="hl str"></span>
<span class="hl str">	perl ./vcf_tab_to_fasta_alignment.pl -i hyp155_n_0.33_mac4_5kb.tab &gt; hyp155_n_0.33_mac4_5kb.fas</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Finally, again `raxml` is run (equivalent to git 13.4).


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 13.7</span>
<span class="hl slc">// run raxml (floridae-rooted)</span>
<span class="hl kwa">process</span> hypo_whg_raxml <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/raxml/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> fas <span class="hl opt">)</span> <span class="hl kwa">from</span> raxml_hypo_genotypes_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;hyp155_n_0.33_mac4_5kb.raxml.support&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> raxml_hypo_whg_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# Infer phylogeny</span>
<span class="hl str">	# Note: number of invariant sites for Felsenstein correction was calculated as number of</span>
<span class="hl str">	# variant sites in alignment (105,043) / genome-wide proportion of variant sites </span>
<span class="hl str">	# (0.05) * genome-wide proportion of invariant sites (0.95)</span>
<span class="hl str"></span>
<span class="hl str">	raxml-NG --all \</span>
<span class="hl str">	  --msa</span> <span class="hl ipl">${fas}</span> <span class="hl str">\</span>
<span class="hl str">	  --model GTR+G+ASC_FELS{1995817} \</span>
<span class="hl str">	  --tree pars{20},rand{20} \</span>
<span class="hl str">	  --bs-trees 100 \</span>
<span class="hl str">	  --threads 24 \</span>
<span class="hl str">	  --worker 8 \</span>
<span class="hl str">	  --seed 123 \</span>
<span class="hl str">	  --prefix hyp155_n_0.33_mac4_5kb</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
