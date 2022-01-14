# Code repository for: *Rapid radiation in a highly diverse marine environment*

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/145688309.svg)](https://zenodo.org/badge/latestdoi/145688309)
<!-- badges: end -->

This repository contains the original bioinformatic analysis behind the paper *Rapid radiation in a highly diverse marine environment* by Hench, Helmkampf, McMillan and Puebla.

It covers all steps from genotyping based on raw sequencing data, over population genetic analysis to the final plotting of the figures used within the publication.

A more extensive documentation of the individual steps can be found  in the `./docs` folder and is also hosted under the accompanying [github page](https://k-hench.github.io/hamlet_radiation/).

There are two more accompanying repositories for this publication:
- The ENA sequencing repository: contains the raw sequencing data (Accession Nr: PRJEB35459)
- The [dryad data repository](https://doi.org/10.5061/dryad.280gb5mmt): contains the genotypes and some intermediate population genetic results

### R setup

There is an additional **R** package needed to run the plotting scripts for the figures ({[GenomicOriginsScripts](https://k-hench.github.io/GenomicOriginsScripts/)}).
This depends on several non-CRAN R-packages, so to be able to install the package successfully, the following packages will also need to be installed:

```r
# installing non-CRAN dependencies
install.packages("remotes")
remotes::install_bioc("rtracklayer")
remotes::install_github("YuLab-SMU/ggtree")
remotes::install_github("k-hench/hypogen")
remotes::install_github("k-hench/hypoimg")
# installing GenomicOriginsScripts
remotes::install_github("k-hench/GenomicOriginsScripts")
```

Once these non-CRAN packages are installed, it should be possible to re-create the used **R** environment using the {[renv](https://rstudio.github.io/renv/)} package.
After opening the RStudio project (`hamlet_radiation.Rproj`), call:

```r
# restoring R environment
install.packages("renv")
renv::restore()
```

### git tags

In case you are looking for specific parts referenced in the *Materials and Methods* section of the study, these refer to the following locations:

- *git  1.x*: [`nf/01_genotyping/genotyping.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/01_genotyping/genotyping.nf)
- *git  2.x*: [`nf/02_genotyping_all_basepairs/genotyping_all_basepairs.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/02_genotyping_all_basepairs/genotyping_all_basepairs.nf)
- *git  3.x*: [`nf/03_analysis_fst_gxp/analysis_fst_gxp.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/03_analysis_fst_gxp/analysis_fst_gxp.nf)
- *git  4.x*: [`nf/04_analysis_dxy/analysis_dxy.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/04_analysis_dxy/analysis_dxy.nf)
- *git  5.x*: [`nf/05_analysis_fasttree_twisst/analysis_fasttree_twisst.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/05_analysis_fasttree_twisst/analysis_fasttree_twisst.nf)
- *git  6.x*: [`nf/06_analysis_recombination/analysis_recombination.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/06_analysis_recombination/analysis_recombination.nf)
- *git  7.x*: [`nf/07_analysis_pca/analysis_pca.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/07_analysis_pca/analysis_pca.nf)
- *git  8.x*: [`nf/08_analysis_msmc/analysis_msmc.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/08_analysis_msmc/analysis_msmc.nf)
- *git  9.x*: [`nf/09_analysis_hybridization/analysis_hybridization.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/09_analysis_hybridization/analysis_hybridization.nf)
- *git 10.x*: [`nf/10_analysis_admixture/analysis_admixture.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/10_analysis_admixture/analysis_admixture.nf)
- *git 11.x*: [`nf/11_analysis_allele_age/analysis_allele_age.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/11_analysis_allele_age/analysis_allele_age.nf)
- *git 12.x*: [`nf/12_analysis_fst_signif/analysis_fst_sign.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/12_analysis_fst_signif/analysis_fst_sign.nf)
- *git 13.x*: [`nf/13_analysis_phylo_whg/analysis_phylo_whg.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/13_analysis_phylo_whg/analysis_phylo_whg.nf)
- *git 14.x*: [`nf/14_analysis_phylo_regions/analysis_phylo_regions.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/14_analysis_phylo_regions/analysis_phylo_regions.nf)
- *git 15.x*: [`nf/15_analysis_pi/analysis_pi.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/15_analysis_pi/analysis_pi.nf)
- *git 16.x*: [`nf/16_analysis_ibd/analysis_ibd.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/16_analysis_ibd/analysis_ibd.nf)
- *git 17.x*: [`nf/17_analysis_dstats/analysis_dstats.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/17_analysis_dstats/analysis_dstats.nf)
- *git 18.x*: [`nf/18_genotyping_all_basepairs_mt/genotyping_all_basepairs_mt.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/18_genotyping_all_basepairs_mt/genotyping_all_basepairs_mt.nf)
- *git 19.x*: [`19_analysis_phylo_serraninae/analysis_phylo_serraninae.nf`](https://github.com/k-hench/hamlet_radiation/blob/master/nf/19_analysis_phylo_serraninae/analysis_phylo_serraninae.nf)
- *git 20.x*: [`sh/create_figures.sh`](https://github.com/k-hench/hamlet_radiation/blob/master/sh/create_figures.sh)

Again, please refer to the [docs](https://k-hench.github.io/hamlet_radiation/) for a more in depth documentation.

---

<p align="center"><img src="logo.svg" alt="logo" width="150"/></p>
