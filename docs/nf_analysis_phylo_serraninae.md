---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 19) Analysis XVI (Serraninae Phylogeny)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/19_analysis_phylo_serraninae
nextflow run analysis_phylo_serraninae.nf -c ../../nextflow.config -resume
```

## Summary

The <span style="color:red;">...</span> are computed within the [**nextflow**](https://www.nextflow.io/) script `analysis_phylo_serraninae.nf` (located under `$BASE_DIR/nf/19_analysis_phylo_serraninae/`).
It takes the <span style="color:red;">...</span> and computes <span style="color:red;">...</span>.
Below is an overview of the steps involved in the analysis.

## Details of `analysis_phylo_serraninae.nf`

### Setup

The nextflow script starts by <span style="color:red;">...</span>

:::kclass

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow

<span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// this pipeline was not actually run using nextflow,</span>
<span class="hl slc">// but managed manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>

<span class="hl slc">// git 19.1</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&#39;../../ressources/serraninae/Serraninae_Rabosky.phy&#39;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> serraninae_phy_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.2</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&#39;../../ressources/serraninae/partitions.txt&#39;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> partitions_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.3</span>
<span class="hl slc">// Obtain gene boundaries and split alignment into individual genes</span>
<span class="hl kwa">process</span> get_gene_boundaries  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;gene_boundaries&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> phy <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> part <span class="hl opt">)</span> <span class="hl kwa">from</span> serraninae_phy_ch.combine<span class="hl opt">(</span> partitions_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*serrR.aln&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> <span class="hl opt">(</span> gene_boundaries_ch<span class="hl opt">,</span> gene_boundaries_ch2 <span class="hl opt">)</span>

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# manually extract Serraninae from Rabosky alignment, final_alignment.phylip (FToL Dryad) &gt; Serraninae_Rabosky.phy</span>
<span class="hl str">	perl \$BASE_DIR/pl/Phylip2Fasta.pl</span> <span class="hl ipl">${phy}</span> <span class="hl str">Serraninae_Rabosky.aln</span>
<span class="hl str"></span>
<span class="hl str">	sed &#39;s/-/N/g&#39; Serraninae_Rabosky.aln &gt; Serraninae_RaboskyN.aln</span>
<span class="hl str"></span>
<span class="hl str">	sed &#39;s/DNA,</span> <span class="hl esc">\\</span><span class="hl str">(.*</span><span class="hl esc">\\</span><span class="hl str">) =</span> <span class="hl esc">\\</span><span class="hl str">(.*</span><span class="hl esc">\\</span><span class="hl str">)-</span><span class="hl esc">\\</span><span class="hl str">(.*</span><span class="hl esc">\\</span><span class="hl str">)/</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">2</span><span class="hl esc">\\</span><span class="hl str">t</span><span class="hl esc">\\</span><span class="hl str">3/g&#39;</span> <span class="hl ipl">${part}</span> <span class="hl str">&gt; partitions.bed   # (partitions.txt from FToL Dryad)</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{ print \$2 }&#39; partitions.bed | tail -n +2 | tr &#39;</span><span class="hl esc">\\</span><span class="hl str">n&#39; &#39;,&#39; | sed &#39;s/.\$//&#39;</span>
<span class="hl str"></span>
<span class="hl str">	# phast from http://compgen.cshl.edu/phast/</span>
<span class="hl str">	\$SFTWR/phast/bin/msa_split Serraninae_RaboskyN.aln -r gene \</span>
<span class="hl str">	  --by-index 980,1757,2292,2974,4115,4955,5681,6614,7244,7988,8822,9869,11402,12142,12952,13663,15117,16341,17265,17910,18633,19728,20715,21522,22350,23115</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{ print \$1, \$2, \$3 }&#39; partitions.bed | xargs -n 3 sh -c &#39;mv gene.\$1-\$2.fa \$0_serrR.aln&#39;</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.4</span>
<span class="hl slc">// Identify genes in reference genome</span>
<span class="hl kwa">process</span> identify_genes  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;identify_genes&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> aln <span class="hl opt">)</span> <span class="hl kwa">from</span> gene_boundaries_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;coord_R24_Hpue_ed.bed&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> gene_coords_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# manually select query sequences from Serraninae_RaboskyN.aln &gt; queries_R24.fas</span>
<span class="hl str"></span>
<span class="hl str">	mkdir blast_db</span>
<span class="hl str">	cp \$BASE_DIR/ressources/HP_genome_unmasked_01.fa blast_db/Hpue_genome_unmasked_01.fas</span>
<span class="hl str">	makeblastdb -in blast_db/Hpue_genome_unmasked_01.fas -dbtype nucl -parse_seqids</span>
<span class="hl str"></span>
<span class="hl str">	blastn -query queries_R24.fas -db blast_db/Hpue_genome_unmasked_01.fas -out blast_R24-Hpue_aln.txt -outfmt 0 -evalue 1e-10</span>
<span class="hl str">	blastn -query queries_R24.fas -db blast_db/Hpue_genome_unmasked_01.fas -out Rblast_R24-Hpue_tab.csv -outfmt 6 -evalue 1e-10</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;OFS = &quot;</span><span class="hl esc">\\</span><span class="hl str">t&quot; { if (\$9 &lt; \$10) print \$2, \$9, \$10, \$1, &quot;+&quot; ; else print \$2, \$10, \$9, \$1, &quot;-&quot; }&#39; blast_R24-Hpue_tab.csv | \</span>
<span class="hl str">		sed &#39;s/</span><span class="hl esc">\\</span><span class="hl str">(.*</span><span class="hl esc">\\</span><span class="hl str">)_</span><span class="hl esc">\\</span><span class="hl str">(.*</span><span class="hl esc">\\</span><span class="hl str">)</span><span class="hl esc">\\</span><span class="hl str">t/</span><span class="hl esc">\\</span><span class="hl str">1</span><span class="hl esc">\\</span><span class="hl str">t/g&#39; &gt; coord_R24_Hpue.csv</span>
<span class="hl str"></span>
<span class="hl str">	# manually combine segmented genes (16s, glyt, tbr1, zic1) in coord_R24_Hpue.csv &gt; coord_R24_Hpue_ed.csv</span>
<span class="hl str"></span>
<span class="hl str">	sed -Ei &#39;s/_[A-Z][a-z]{3}//g&#39; coord_R24_Hpue_ed.bed</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.5</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&#39;../../1_genotyping/2_raw_vcfs/all_sites.unplaced.vcf.gz&#39;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_unlplaced_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.6</span>
<span class="hl slc">// Prepare contigs </span>
<span class="hl kwa">process</span> prepare_contigs  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;prepare_contigs&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span> <span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_unlplaced_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_sites.Contig*.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.vcf.gz.tbi&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> contigs_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	zcat &lt;</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| grep -e &#39;#&#39; -e &#39;Contig11544&#39; | bgzip &gt; all_sites.Contig11544.vcf.gz</span>
<span class="hl str">	zcat &lt;</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| grep -e &#39;#&#39; -e &#39;Contig11607&#39; | bgzip &gt; all_sites.Contig11607.vcf.gz</span>
<span class="hl str">	zcat &lt;</span> <span class="hl ipl">${vcf}</span> <span class="hl str">| grep -e &#39;#&#39; -e &#39;Contig11888&#39; | bgzip &gt; all_sites.Contig11888.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	tabix -p vcf *.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.7</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&#39;../../ressources/serraninae/samples_to_include.ids&#39;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> samples_ch <span class="hl opt">}</span>
</code>
</pre>
</div>



<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.8</span>
<span class="hl slc">// Extract genotypes and convert to Fasta</span>
<span class="hl kwa">process</span> extract_genotypes  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;extract_genotypes&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> crds <span class="hl opt">),</span> <span class="hl kwc">file</span> <span class="hl opt">(</span> samples <span class="hl opt">)</span> <span class="hl kwa">from</span> gene_coords_ch.combine<span class="hl opt">(</span> samples_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> genotypes_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# select representative samples based on highest coverage (except outliers) and regional diversity &gt; samples_to_include.ids</span>
<span class="hl str"></span>
<span class="hl str">	awk &#39;{ print $1, $2, $3, $4 }&#39;</span> <span class="hl ipl">${crds}</span> <span class="hl str">| \</span>
<span class="hl str">		xargs -n 4 sh -c &#39;vcftools --gzvcf \$BASE_DIR/2_raw_vcfs/all_sites.&quot;\$0&quot;.vcf.gz --keep</span> <span class="hl ipl">${samples}</span> <span class="hl str">--chr &quot;\$0&quot; --from-bp &quot;\$1&quot; --to-bp &quot;\$2&quot; --recode --stdout | grep -v &quot;##&quot; &gt; &quot;\$3&quot;_hypS.vcf&#39;</span>
<span class="hl str"></span>
<span class="hl str">	for FILE in ./*.vcf</span>
<span class="hl str">		do</span>
<span class="hl str">		vcf-to-tab &lt; \$FILE &gt; \${FILE%.vcf}.tab</span>
<span class="hl str">		perl \$SFTWR/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i \${FILE%.vcf}.tab &gt; \${FILE%.vcf}.fas</span>
<span class="hl str">		rm \${FILE%.vcf}.tab*</span>
<span class="hl str">	done</span>
<span class="hl str"></span>
<span class="hl str">	rev=`awk &#39;\$5 == &quot;-&quot; { print \$4 }&#39;</span> <span class="hl ipl">${crds}</span><span class="hl str">`</span>
<span class="hl str"></span>
<span class="hl str">	for GENE in \$rev</span>
<span class="hl str">		do</span>
<span class="hl str">		seqtk seq -r \</span><span class="hl ipl">${GENE}</span><span class="hl str">_hypS.fas &gt; \</span><span class="hl ipl">${GENE}</span><span class="hl str">_hypS_rev.fas</span>
<span class="hl str">		mv \</span><span class="hl ipl">${GENE}</span><span class="hl str">_hypS.fas \</span><span class="hl ipl">${GENE}</span><span class="hl str">_hypS.fas_org</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>

<span class="hl slc">// git 19.9</span>
<span class="hl kwa">Channel</span>
	.fromPath<span class="hl opt">(</span><span class="hl str">&#39;../../ressources/serraninae/taxa_to_exclude.ids&#39;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> taxa_ch <span class="hl opt">}</span>

<span class="hl slc">// git 19.10</span>
<span class="hl slc">// Combine Rabosky and present data</span>
<span class="hl kwa">process</span> combine_data  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;combine_data&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> aln <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> taxa <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> fa <span class="hl opt">)</span> <span class="hl kwa">from</span> gene_boundaries_ch2.combine<span class="hl opt">(</span> taxa_ch <span class="hl opt">)</span>.combine<span class="hl opt">(</span> genotypes_ch <span class="hl opt">)</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_new.fas&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> data_combined_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	# manually prepared list of taxa to exclude from Rabosky&#39;s alignments &gt; taxa_to_exclude.ids</span>
<span class="hl str"></span>
<span class="hl str">	for FILE in *.aln</span>
<span class="hl str">		do</span>
<span class="hl str">		sed -e &#39;s/N//g&#39; -e &#39;s/&gt; /&gt;/g&#39; -e &#39;/^\$/d&#39; \$FILE | </span>
<span class="hl str">		awk &#39;BEGIN { while ((getline &lt; &quot;</span><span class="hl ipl">${taxa}</span><span class="hl str">&quot;) &gt;0) l[&quot;&gt;&quot;\$1]=1 } /^&gt;/ { f=!l[\$1] }f&#39; &gt; \${FILE%.aln}_pre.fa</span>
<span class="hl str">	done</span>
<span class="hl str"></span>
<span class="hl str">	for FILE in *.fas</span>
<span class="hl str">		do</span>
<span class="hl str">		NAME=\${FILE%_hypS.fas}</span>
<span class="hl str">		cat \$FILE \</span><span class="hl ipl">${NAME}</span><span class="hl str">_serrR_pre.fa &gt; \</span><span class="hl ipl">${NAME}</span><span class="hl str">_new.fas</span>
<span class="hl str">	done</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// ----------------------- DISCLAIMER ----------------------</span>
<span class="hl slc">// Between git 19.11 and git 19.12 the alignments </span>
<span class="hl slc">// were additionally currated manually</span>
<span class="hl slc">// ---------------------------------------------------------</span>
</code>
</pre>
</div>


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.11</span>
<span class="hl slc">// Align sequences and remove ambiguously aligned positions</span>
<span class="hl kwa">process</span> align_sequences  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;align_sequences&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> fa <span class="hl opt">)</span> <span class="hl kwa">from</span> data_combined_ch

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.aln&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> align_sequences_ch

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	for FILE in ./*_new.fas</span>
<span class="hl str">		do </span>
<span class="hl str">		mafft --maxiterate 1000 --globalpair --adjustdirectionaccurately \</span><span class="hl ipl">${FILE}</span> <span class="hl str">&gt; \${FILE%.fas}.aln</span>
<span class="hl str">	done</span>
<span class="hl str"></span>
<span class="hl str">	Gblocks 3_alignments/12s_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100</span>
<span class="hl str">	Gblocks 3_alignments/16s_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100</span>
<span class="hl str">	Gblocks 3_alignments/coi_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100</span>
<span class="hl str">	Gblocks 3_alignments/cytb_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100</span>
<span class="hl str"></span>
<span class="hl str">	# finalize with manual alignment editing step:</span>
<span class="hl str">	# 12S:  removed 45±3 bp (whole gap)           &gt; 12s_new_ed.aln</span>
<span class="hl str">	# 16S:  removed positions 194, 252, 315       &gt; 16s_new_ed.aln</span>
<span class="hl str">	# coi:  removed entire gene in torpan, tabhon &gt; coi_new_ed.aln</span>
<span class="hl str">	# cytb: removed entire gene in torpan, tabhon &gt; cytb_new_ed.aln</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 19.12</span>
<span class="hl slc">// Phylogenetic inference (IQ-TREE)</span>
<span class="hl kwa">process</span> phylogenetic_inference  <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;phylogenetic_inference&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../2_analysis/fotl/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span> 

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> aln <span class="hl opt">)</span> <span class="hl kwa">from</span> align_sequences_ch
	
	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;concat_R24ed.treefile&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> phylo_out

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	mkdir alignments</span>
<span class="hl str">	cp *.aln alignments/</span>
<span class="hl str">	iqtree2 -p ./alignments --prefix concat_R24ed -B 1000 -T AUTO</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
:::

---
