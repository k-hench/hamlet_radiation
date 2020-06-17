---
title: "Script repository"
subtitle: "(Hench *et al.* supplement)"
author: "Kosmas Hench"
date: "2020-06-17"
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: k-hench/bookdown
description: "Scripts used to produce Figures and Supplementary Figures of 'Ancestral variation, hybridization and modularity fuel a marine radiation' by Hench, McMillan an Puebla"
---

# Intro



This repository contains the complete workflow used in the paper "*Ancestral variation, hybridization and modularity fuel a marine radiation*".
The individual chapters of this documentation follow the separate main steps of the workflow, which each refer to an individual prefix in the _git x.x_ references of the papers method section.
The individual steps partly depend on each other - especially git 1 - git 3 should be executed in order and before the other steps.

<!--<div style="max-width:800px; margin:auto;">
<img src="index_files/figure-html/unnamed-chunk-1-1.png" width="864" />
</div> --->

## Analysis

A documentation of the data preparation and the data analysis (git 1.x - 9.x) can be found at:

- git 1.x: [Genotyping](genotyping-i-snps-only.html)
- git 2.x: [Genotyping all base pairs](genotyping-ii-all-callable-sites.html)
- git 3.x: [Analysis (<i>F<sub>ST</sub></i> & GxP)](analysis-i-fst-gxp.html)
- git 4.x: [Analysis (<i>d<sub>XY</sub></i> & $\pi$)](analysis-ii-dxy-pi.html)
- git 5.x: [Analysis (phylogeny & topolgy weighting)](analysis-iii-phylogeny-topology-weighting.html)
- git 6.x: [Analysis ($\rho$)](analysis-iv-rho.html)
- git 7.x: [Analysis (<i>H<sub>O</sub></i> )](analysis-v-ho.html)
- git 8.x: [Analysis (demographic history)](analysis-vi-demographic-history.html)
- git 10.x: [Analysis (Stankowski *et al.* 2019)](analysis-vii-monkeyflowers.html)

## Prerequesites

All scripts assume two variables to be set within the bash environment:

  - `$BASE_DIR` is assumed to point to the base folder of this repository
  - `$SFTWR` is a folder that contains all the software dependencies that are used within the scripts

Furthermore, external dependencies need to be downloaded and deployed at the expected places (s. README.md at the `ressources` folder).

## Figures

The creation of the figures is bundled in a single script (git 11) which can be executed once all `nextflow` scripts have successfully run.

```sh
cd $BASE_DIR
bash sh/create_figures.sh
```

This is basically just a collection that will run all scripts located under `$BASE_DIR/R/fig`.
Under this location, you will find one `R` script per figure (and suppl. figure).
So if you are only interested in a single figure - that is the place to start looking.

Furthermore, a more detailed documentation exists for all the figure scripts used for the manuscript:

[F1](figure-1.html), [F2](figure-2.html), [F3](figure-3.html) & [F4](figure-4.html)

as well as for all the supplementary figures:

[SF1](), [SF2](), [SF3](), [SF4](), [SF5](), [SF6](), [SF7](), [SF8]() and [SF9]().

