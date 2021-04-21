---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---



# (git 14) Analysis XII  (Outlier Region Phylogenies)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/14_analysis_phylo_regions
nextflow run analysis_phylo_regions.nf
```

## Summary

The phylogenies specific to particular differentiation outlier regions are reconstructed within the [**nextflow**](https://www.nextflow.io/) script `analysis_phylo_regions.nf` (located under `$BASE_DIR/nf/14_analysis_phylo_regions/`).
This includes both the *sample-level* as well as the *population-level* phylogenies.

## Details of `analysis_phylo_regions.nf`

> This part of the analysis was actually manged manually and not via `nextflow`. 
> We still report the analysis as a `.nf` script as we believe this is a cleaner and more concise report of the conducted analysis.

### Setup

The nextflow script starts by opening the two specific linkage groups of the *all_bp* genotype data set and binding it to differentiation outlier IDs as well as to a reference table containing the genomic coordinates of all differentiation outlier regions.

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow

<span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// this pipeline was not actually run using nexflow,</span>
<span class="hl slc">// but managed manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>

<span class="hl slc">// Region-specific phylogenies</span>
<span class="hl slc">// ---------------------------</span>

<span class="hl slc">// git 14.1</span>
<span class="hl slc">// bundle allBP files and outlier table</span>
<span class="hl kwa">Channel</span>
  .fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG04.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
  .concat<span class="hl opt">(</span><span class="hl kwa">Channel</span>.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">))</span>
  .concat<span class="hl opt">(</span><span class="hl kwa">Channel</span>.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/byLG/filterd.allBP.LG12.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">))</span>
  .merge<span class="hl opt">(</span><span class="hl kwa">Channel</span>.from<span class="hl opt">(</span><span class="hl str">&quot;LG04_1&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;LG12_3&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;LG12_4&quot;</span><span class="hl opt">))</span>
  .combine<span class="hl opt">(</span><span class="hl kwa">Channel</span>.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/focal_outlier.tsv&quot;</span><span class="hl opt">))</span>
  .set<span class="hl opt">{</span> vcf_lg_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Then, two different sample lists (both excluding hybrid samples, one with and one without Serranus outgroup samples) are loaded and bound to a sample mode identifier.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.2</span>
<span class="hl slc">// toggle sample modes (with / without Serranus outgroup)</span>
<span class="hl kwa">Channel</span>.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/samples_155.txt&quot;</span><span class="hl opt">)</span>
  .concat<span class="hl opt">(</span><span class="hl kwa">Channel</span>.fromPath<span class="hl opt">(</span><span class="hl str">&quot;../../ressources/samples_hybrids.txt&quot;</span><span class="hl opt">))</span>
  .merge<span class="hl opt">(</span><span class="hl kwa">Channel</span>.from<span class="hl opt">(</span><span class="hl str">&quot;155&quot;</span><span class="hl opt">,</span> <span class="hl str">&quot;hyS&quot;</span><span class="hl opt">))</span>
  .set<span class="hl opt">{</span> sample_mode_ch <span class="hl opt">}</span>
</code>
</pre>
</div>

Next, for each outlier region, the genotype data is subset to the respective outlier region.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.3</span>
<span class="hl slc">// subset genotypes to outlier region</span>
<span class="hl kwa">process</span> extract_regions <span class="hl opt">{</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> vcfIdx <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> outlier_file <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> sample_file <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_lg_ch.combine<span class="hl opt">(</span> sample_mode_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.vcf&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> vcf_raxml_ch<span class="hl opt">,</span> vcf_pomo_ch <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# Extract regions of interest from genotype data (allBP),</span>
<span class="hl str">	# remove hybrid / Serranus samples and indels; simplify headers</span>
<span class="hl str"></span>
<span class="hl str">	head -n 1</span> <span class="hl ipl">${outlier_file}</span> <span class="hl str">| cut -f 1-3 &gt; outlier.bed</span>
<span class="hl str">	grep</span> <span class="hl ipl">${outlierId} ${outlier_file}</span> <span class="hl str">| cut -f 1-3 &gt;&gt; outlier.bed</span>
<span class="hl str"></span>
<span class="hl str">	OUT_ALT=\$(echo</span> <span class="hl ipl">${outlierId}</span> <span class="hl str">| tr &#39;[:upper:]&#39; &#39;[:lower:]&#39; | sed &#39;s/_/./&#39;)</span>
<span class="hl str"></span>
<span class="hl str">	vcftools --gzvcf \</span>
<span class="hl str"></span>	  <span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">	  --bed outlier.bed \</span>
<span class="hl str">	  --remove-indels \</span>
<span class="hl str">	  --remove</span> <span class="hl ipl">${sample_file}</span> <span class="hl str">\</span>
<span class="hl str">	  --recode \</span>
<span class="hl str">	  --stdout | \</span>
<span class="hl str">	  grep -v &#39;##&#39; &gt; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.vcf</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then, for the *population-level* phylogenies, the genotypes are first converted to `fasta` format and then to a allele frequency format (`.cf`).
At that point, `iqtree2` is run to create the *population-level* phylogenies.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.4</span>
<span class="hl slc">// run iqtree under pomo model</span>
<span class="hl kwa">process</span> run_pomo <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/revPoMo/outlier_regions/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_raxml_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_pop.cf.treefile&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> pomo_results_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	OUT_ALT=\$(echo</span> <span class="hl ipl">${outlierId}</span> <span class="hl str">| tr &#39;[:upper:]&#39; &#39;[:lower:]&#39; | sed &#39;s/_/./&#39;)</span>
<span class="hl str"></span>
<span class="hl str">	# Convert to fasta format (Python scripts available at https://github.com/simonhmartin/genomics_general), picked up from 6.1.1 output</span>
<span class="hl str">	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py -i</span> <span class="hl ipl">${vcf}</span> <span class="hl str">&gt; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.geno</span>
<span class="hl str"></span>
<span class="hl str">	python \$SFTWR/genomics_general/genoToSeq.py \</span>
<span class="hl str">		-g \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.geno \</span>
<span class="hl str">		-s \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.fas \</span>
<span class="hl str">		-f fasta \</span>
<span class="hl str">		--splitPhased</span>
<span class="hl str">	</span>
<span class="hl str">	# Reformat sample ids to provide population prefixes for cflib</span>
<span class="hl str">	sed -e &#39;s/-/_/g&#39; -e &#39;s/&gt;\(.*\)\([a-z]\{6\}\)_\([AB]\)/&gt;\2-\1_\3/g&#39; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">.fas &gt; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">_p.fas</span>
<span class="hl str"></span>
<span class="hl str">	# Convert to allele frequency format (cflib library available at https://github.com/pomo-dev/cflib)</span>
<span class="hl str">	\$SFTWR/cflib/FastaToCounts.py \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">_p.fas \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">_pop.cf</span>
<span class="hl str"></span>
<span class="hl str">	# IQTREE analysis under PoMo model</span>
<span class="hl str">	iqtree2 \</span>
<span class="hl str">		-nt 16 \</span>
<span class="hl str">		-s \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">_pop.cf \</span>
<span class="hl str">		-m HKY+F+P+N9+G4 \</span>
<span class="hl str">		-b 100</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

For the *sample-level* phylogenies, the genotypes are also converted to `fasta` format.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.5</span>
<span class="hl slc">// convert genotypes to fasta for raxml</span>
<span class="hl kwa">process</span> conversion_raxml <span class="hl opt">{</span>
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_pomo_ch
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*N.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> outlier_regions_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	OUT_ALT=\$(echo</span> <span class="hl ipl">${outlierId}</span> <span class="hl str">| tr &#39;[:upper:]&#39; &#39;[:lower:]&#39; | sed &#39;s/_/./&#39;)</span>
<span class="hl str"></span>
<span class="hl str">	# Replace unknown character states and asterisks (deletions as encoded by GATK) with &quot;N&quot;</span>
<span class="hl str">	vcf-to-tab &lt;</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| sed -e &#39;s/</span><span class="hl esc">\\</span><span class="hl str">.</span><span class="hl esc">\\</span><span class="hl str">/</span><span class="hl esc">\\</span><span class="hl str">./N</span><span class="hl esc">\\</span><span class="hl str">/N/g&#39; -e &#39;s/[ACGTN</span><span class="hl esc">\\</span><span class="hl str">*]</span><span class="hl esc">\\</span><span class="hl str">/</span><span class="hl esc">\\</span><span class="hl str">*/N</span><span class="hl esc">\\</span><span class="hl str">/N/g&#39; &gt; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">N.tab</span>
<span class="hl str"></span>
<span class="hl str">	# Convert to fasta format (Perl script available at https://github.com/JinfengChen/vcf-tab-to-fasta)</span>
<span class="hl str">	wget https://raw.githubusercontent.com/JinfengChen/vcf-tab-to-fasta/master/vcf_tab_to_fasta_alignment.pl</span>
<span class="hl str">	perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">N.tab &gt; \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">N.fas</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

Then, `raxml` is run directly on the `fasta` files.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 14.6</span>
<span class="hl slc">// run raxml</span>
<span class="hl kwa">process</span> run_raxml <span class="hl opt">{</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/raxml/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 
	
	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> outlierId <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> sample_mode <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> fas <span class="hl opt">)</span> <span class="hl kwa">from</span> outlier_regions_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.raxml.support&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> outlier_results_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	OUT_ALT=\$(echo</span> <span class="hl ipl">${outlierId}</span> <span class="hl str">| tr &#39;[:upper:]&#39; &#39;[:lower:]&#39; | sed &#39;s/_/./&#39;)</span>
<span class="hl str"></span>
<span class="hl str">	# Reconstruct phylogenies</span>
<span class="hl str">	raxml-NG --all \</span>
<span class="hl str">		--msa</span> <span class="hl ipl">${fas}</span> <span class="hl str">\</span>
<span class="hl str">		--model GTR+G \</span>
<span class="hl str">		--tree pars{10},rand{10} \</span>
<span class="hl str">		--bs-trees 100 \</span>
<span class="hl str">		--threads 24 \</span>
<span class="hl str">		--worker 8 \</span>
<span class="hl str">		--seed 123 \</span>
<span class="hl str">		--prefix \</span><span class="hl ipl">${OUT_ALT}</span><span class="hl str">_</span><span class="hl ipl">${sample_mode}</span><span class="hl str">N</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
