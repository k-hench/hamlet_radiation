# Code repository for: *Ancestral variation, hybridization and modularity fuel a marine radiation*

This repository contains the original bioinformatic analysis behind the paper *Ancestral variation, hybridization and modularity fuel a marine radiation* Hench, Helmkampf, McMillan and Puebla.

It covers all steps from genotyping based on raw sequencing data, over population genetic analysis to the final plotting of the figures and compilation of tables used within the publication.

A more extensive documentation of the individual steps can be found  in the `./docs` folder and is also hosted under the accompanying [github page](https://k-hench.github.io/hamlet_radiation/).

There are two more accompanying repositories for this publication:
- The ENA sequencing repository: contains the raw sequencing data (Accession Nr: PRJEB35459)
- The [dryad data repository](https://doi.org/10.5061/dryad.280gb5mmt): contains the genotypes and some intermediate population genetic results

There is an additional [**R** package](https://k-hench.github.io/GenomicOriginsScripts/) needed to run the plotting scripts for the figures.
GenomicOriginsScripts depends on several non-CRAN R-packages.
To be able to install the package successfully, the following packages will also need to be installed:

```r
# installing dependencies
install.packages("remotes")
remotes::install_bioc("rtracklayer")
remotes::install_github("k-hench/hypogen")
remotes::install_github("k-hench/hypoimg")
# installing GenomicOriginsScripts
remotes::install_github("k-hench/GenomicOriginsScripts")
```

In case you are looking for specific parts referenced in the *Materials and Methods* section of the study, these refer to the following locations:

- *git  1*: [`nf/01_genotyping/genotyping.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/01_genotyping/genotyping.nf)
- *git  2*: [`nf/02_genotyping_all_basepairs/genotyping_all_basepairs.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/02_genotyping_all_basepairs/genotyping_all_basepairs.nf)
- *git  3*: [`nf/03_analysis_fst_gxp/analysis_fst_gxp.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/03_analysis_fst_gxp/analysis_fst_gxp.nf)
- *git  4*: [`nf/04_analysis_dxy/analysis_dxy.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/04_analysis_dxy/analysis_dxy.nf)
- *git  5*: [`nf/05_analysis_fasttree_twisst/analysis_fasttree_twisst.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/05_analysis_fasttree_twisst/analysis_fasttree_twisst.nf)
- *git  6*: [`nf/06_analysis_recombination/analysis_recombination.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/06_analysis_recombination/analysis_recombination.nf)
- *git  7*: [`nf/07_analysis_pca/analysis_pca.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/07_analysis_pca/analysis_pca.nf)
- *git  8*: [`nf/08_analysis_msmc/analysis_msmc.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/08_analysis_msmc/analysis_msmc.nf)
- *git  9*: [`nf/09_analysis_hybridization/analysis_hybridization.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/09_analysis_hybridization/analysis_hybridization.nf)
- *git 10*: [`nf/10_analysis_admixture/analysis_admixture.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/10_analysis_admixture/analysis_admixture.nf)
- *git 12*: [`nf/11_analysis_allele_age/analysis_allele_age.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/11_analysis_allele_age/analysis_allele_age.nf)
- *git 13*: [`nf/12_analysis_fst_signif/analysis_fst_sign.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/12_analysis_fst_signif/analysis_fst_sign.nf)
- *git 14*: [`nf/14_analysis_phylo_regions/analysis_phylo_regions.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/14_analysis_phylo_regions/analysis_phylo_regions.nf)
- *git 15*: [`sh/create_figures.sh`](https://github.com/k-hench/hamlet_radiation/blob/master/sh/create_figures.sh)

Again, please refer to the [docs](https://k-hench.github.io/hamlet_radiation/) for a more in depth documentation.

---

<p align="center"><img src="logo.svg" alt="logo" width="150"/></p>
