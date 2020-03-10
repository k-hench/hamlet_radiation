# The **R** folder

This folder contains analysis run in **R**.

In priciple, the analysis should be run from within the **R**-project (`*.Rproj`) contained in the repository root folder or they are run from within the **nextflow** pipelines in the root folder.

The scripts to create all figures are contained in the sub-folder `./fig` and are executable R-scripts that should be run from the root folder of the project, eg:

```sh
Rscript --vanilla R/fig/plot_XX.R input_parameter1 input_parameter2
```

Please refer to the [documentation](https://k-hench.github.io/chapter2/) for some more in depth explanations.

Some of the packages need *non-standard* installation, eg.:

- **hypogen**: `remotes::install_github("k-hench/hypogen")`
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
remotes::install_github("hadley/tidyverse")
remotes::install_github("thomasp85/ggforce")
remotes::install_github("slowkow/ggrepel")
remotes::install_github("andland/logisticPCA")
remotes::install_github("sjp/grImport2")
remotes::install_github("k-hench/hpogen")
remotes::install_github("k-hench/hypoimg")
remotes::install_github("k-hench/geomfactory")
remotes::install_github("k-hench/GenomicOriginsScripts")
BiocManager::install("rtracklayer", version = "3.8")
BiocManager::install("SNPRelate", version = "3.8")
```
---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>