# The **R** folder

This folder contains analysis run in **R**.

In priciple, the analysis should be run from within the **R**-project (`*.Rproj`) contained in the repository root folder or they are run from within the **nextflow** pipelines in the root folder.

Some of the packages need *non-standard* installation:

- **ggforce** `devtools::install_github("thomasp85/ggforce")`
- **scico**: `devtools::install_github("thomasp85/scico")`
- **grConvert**: `devtools::install_github('sjp/grConvert')`
- **hypogen**: `devtools::install_github("k-hench/hypogen")`
- **hypoimg**: `devtools::install_github("k-hench/hypoimg")`
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

---

<center><img src="../logo.svg" alt="logo" width="150"/></center>