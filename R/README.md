# The **R** folder

This folder contains analysis run in **R**.

In priciple, the analysis should be run from within the **R**-project (`*.Rproj`) contained in the repository root folder. 

The **R** packages of the project are managed using [**packrat**](https://rstudio.github.io/packrat/).

Some of the packages need *non-standard* installatino though:

- **pcadapt**: `devtools::install_github("bcm-uga/pcadapt")`
- **ggforce** `devtools::install_github("thomasp85/ggforce")`
- **scico**: `devtools::install_github("thomasp85/scico")`
- **grConvert**: `devtools::install_github('sjp/grConvert')`
- **hypogen**: `devtools::install_github("k-hench/hypogen")`
- **hypoimg**: `devtools::install_github("k-hench/hypoimg")`
- **rtracklayer**:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
```

---

<center><img src="../logo.svg" alt="logo" width="150"/></center>