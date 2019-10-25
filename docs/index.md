---
title: "Script repository"
subtitle: "(Hench *et al.* supplement)"
author: "Kosmas Hench"
date: "2019-10-25"
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: k-hench/bookdown
description: "Scripts used to produce Figures and Supplementary Figures of 'The genomic origins of a marine radiation' by Hench, McMillan an Puebla"
---

# Intro

This repository contains a collection of scripts used in the paper "*The genomic origins of a marine radiation*".

## Analysis

A documentation of the data preparation and the data analysis can be found at:

[Genotyping](genotyping.html)

## Figures

A more detailed documentation exists for all the figures of the manuscript:

[F1](figure-1.html), [F2](figure-2.html) & [F3](figure-3.html) 

as well as for all the supplementary figures:

All scripts assume two variables to be set within the bash environment:

  - `$BASE_DIR` is assumed to point to the base folder of this repository
  - `$SFTWR` is a folder that contains all the software dependencies that are used within the scripts

The dependencies need to be downloaded and installed separately.


