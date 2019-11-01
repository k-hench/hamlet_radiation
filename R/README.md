# The **R** folder

This folder contains analysis run in **R**.

In priciple, the analysis should be run from within the **R**-project (`*.Rproj`) contained in the repository root folder or they are run from within the **nextflow** pipelines in the root folder.

Some of the packages need *non-standard* installation, eg.:

- **hypogen**: `devtools::install_github("k-hench/hypogen")`
- **rtracklayer** & **SNPRelate**: (depending on **R** version)
```R
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
biocLite("SNPRelate")
```
or
```
BiocManager::install("rtracklayer", version = "3.8")
BiocManager::install("SNPRelate", version = "3.8")
```

Here is what I ran to set up **R3.5.2**:

```R
install.packages("devtools")
install.packages("BiocManager")
devtools::install_github("hadley/tidyverse")
devtools::install_github("thomasp85/ggforce")
devtools::install_github("slowkow/ggrepel")
devtools::install_github("andland/logisticPCA")
devtools::install_github("sjp/grImport2")
devtools::install_github("k-hench/hpogen")
devtools::install_github("k-hench/hypoimg")
devtools::install_github("k-hench/geomfactory")
devtools::install_github("k-hench/GenomicOriginsScripts")
BiocManager::install("rtracklayer", version = "3.8")
BiocManager::install("SNPRelate", version = "3.8")
```
---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>