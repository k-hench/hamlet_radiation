# The **chapter2** repository

This project contains the analysis of differentiation, divergence and  the phylogenetic relationships between hamlet species (*Hypoplectrus* spp) from Belize, Honduras and Panama.

Steps involved include:

- Genotyping from raw sequencing reads (`*.fq.gz` ==> `phased_mac2.vcf.gz`) using [**GATK4**](https://software.broadinstitute.org/gatk/) and phasing using [**shapeit2**](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- Analysis of differentiation and divergence (*F<sub>ST</sub>*, *D<sub>xy</sub>*, G x P) using [**VCFtools**](https://vcftools.github.io/examples.html), **popgenWindows.py** (by [Simon Martin](https://github.com/simonhmartin/genomics_general)) and [**GEMMA**](https://github.com/genetics-statistics/GEMMA).
- Phylogenetic analysis using [**PhyML**](http://www.atgc-montpellier.fr/index.php?type=bn)/[**Fasttree2**](http://www.microbesonline.org/fasttree/) in combination with [**twisst**](https://github.com/simonhmartin/twisst).
- Analysis of demographic history using [**msmc2**](https://github.com/stschiff/msmc2).

The whole process is monitored using [**Nextflow**](https://www.nextflow.io/index.html) and distributed over several pipelines.
The Nexflow pipelines are stored in individual sub folders within the nf folder:

- genotyping (including variant calling and phasing) is done in `.nf/genotyping/genotyping.nf` and `./nf/genotyping_all_basepairs/genotyping_all_basepairs.nf`.
- population genetic analysis happens within the `./nf/analysis_*/analysis*.nf` pipelines.

The commands to run the pipelines (`nf_run_gatk`) can be found in `sh/nextflow_alias.sh` which are always executed from within the respective analysis folder (eg `nf_run_dxy` is called from `nf/analysis_dxy`).

Additional information about the input data is located in the `./metadata` folder,
a summary of the nextflow processes can be found in the `./docs` folder and is also hosted under the accompanying [github page](https://k-hench.github.io/chapter2/).

There is an additional [**R** package](https://k-hench.github.io/GenomicOriginsScripts/) needed to run the plotting scripts for the figures.
It can be installed using:

```r
remotes::install_github("k-hench/GenomicOriginsScripts")
```

---

<p align="center"><img src="logo.svg" alt="logo" width="150"/></p>
