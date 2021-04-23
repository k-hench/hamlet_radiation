---
title: "Script repository"
subtitle: "(Hench *et al.* supplement)"
author: "Kosmas Hench and Martin Helmkampf"
date: "2021-04-23"
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: k-hench/bookdown
description: "Scripts used to produce Figures and Supplementary Figures of 'Ancestral variation, hybridization and modularity fuel a marine radiation' by Hench, Helmkampf, McMillan an Puebla"
---

# Intro



[This repository](https://github.com/k-hench/chapter2) contains the complete workflow used in the paper "*Ancestral variation, hybridization and modularity fuel a marine radiation*".
The individual chapters of this documentation follow the separate main steps of the workflow.
Each of the chapters thus refers to an individual prefix in the _git x.x_ references of the papers method section.
The individual steps partly depend on each other - especially git 1 - git 3 should be executed in order and before the other steps.

## Analysis

A documentation of the data preparation and the data analysis (git 1.x - 14.x) can be found at:

- git 1.x: [Genotyping](git-1-genotyping-i-snps-only.html)
- git 2.x: [Genotyping all base pairs](git-2-genotyping-ii-all-callable-sites.html)
- git 3.x: [Analysis (<i>F<sub>ST</sub></i> & GxP)](git-3-analysis-i-fst-gxp.html)
- git 4.x: [Analysis (<i>d<sub>XY</sub></i> & $\pi$)](git-4-analysis-ii-dxy-pi.html)
- git 5.x: [Analysis (topolgy weighting)](git-5-analysis-iii-topology-weighting.html)
- git 6.x: [Analysis ($\rho$)](git-6-analysis-iv-rho.html)
- git 7.x: [Analysis (PCA)](git-7-analysis-v-principal-component-analysis.html)
- git 8.x: [Analysis (demographic history)](git-8-analysis-vi-demographic-history.html)
- git 9.x: [Analysis (hybridization)](git-9-analysis-vii-hybridization.html)
- git 10.x: [Analysis (admixture)](git-10-analysis-viii-admixture.html)
- git 11.x: [Analysis (allele age)](git-11-analysis-ix-allele-age.html)
- git 12.x: [Analysis (<i>F<sub>ST</sub></i> permutation)](git-12-analysis-x-fst-permutation-test.html)
- git 13.x: [Analysis (whg phylogeny)](git-13-analysis-xi-whole-genome-phylogenies.html)
- git 14.x: [Analysis (outlier region phylogeny)](git-14-analysis-xii-outlier-region-phylogenies.html)

## Prerequesites

All scripts assume two variables to be set within the bash environment:

  - `$BASE_DIR` is assumed to point to the base folder of this repository
  - `$SFTWR` is a folder that contains all the software dependencies that are used within the scripts

The analysis is controlled using the workflow manager [`nextflow`](https://www.nextflow.io/) and uses slightly different configurations across the individual pipelines. The exact commands used to execute the analysis during the development of the publication are stored within the aliases set within `sh/nextflow_alias.sh`.

Furthermore, external dependencies need to be downloaded and deployed at the expected places (s. README.md at the `ressources` folder).

## Figures

The creation of the figures is bundled in a single script (git 15) which can be executed once all `nextflow` scripts have successfully run.

```sh
cd $BASE_DIR
bash sh/create_figures.sh
```

This is basically just a wrapper script that will run all scripts located under `$BASE_DIR/R/fig`.
Under this location, you will find one `R` script per figure (and suppl. figure).
So if you are only interested in a single figure - that is the place to start looking.

Furthermore, a more detailed documentation exists for all the figure scripts used for the manuscript:

[F1](figure-1.html), [F2](figure-2.html), [F3](figure-3.html), [F4](figure-4.html) [F5](figure-5.html) and [F6](figure-6.html)

as well as for all the supplementary figures:

[SF1](supplementary-figure-1.html), [SF2](supplementary-figure-2.html), [SF3](supplementary-figure-3.html),
[SF4](supplementary-figure-4.html), [SF5](supplementary-figure-5.html), [SF6](supplementary-figure-6.html),
[SF7](supplementary-figure-7.html), [SF8](supplementary-figure-8.html), [SF9](supplementary-figure-9.html),
[SF10](supplementary-figure-10.html), [SF11](supplementary-figure-11.html), [SF12](supplementary-figure-12.html),
[SF13](supplementary-figure-13.html), [SF14](supplementary-figure-14.html), [SF15](supplementary-figure-15.html)
and [SF16](supplementary-figure-16.html).

## R setup

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

Apart from the specific **R** packages that can be retrieved via {renv} from the `renv.lock` file, the used **R** setup at the at time of compilation is as follows:


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.2 LTS
## 
## Matrix products: default
## BLAS:   /usr/local/lib/R/lib/libRblas.so
## LAPACK: /usr/local/lib/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] compiler_4.0.3    magrittr_2.0.1    bookdown_0.19     htmltools_0.5.1.1
##  [5] tools_4.0.3       yaml_2.2.1        stringi_1.5.3     rmarkdown_2.7.6  
##  [9] knitr_1.31        stringr_1.4.0     digest_0.6.27     xfun_0.22        
## [13] rlang_0.4.10      evaluate_0.14
```

---
