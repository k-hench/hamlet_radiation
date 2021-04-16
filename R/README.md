# The **R** folder

This folder contains analysis run in **R**.

In priciple, the analysis should be run from within the **R**-project (`*.Rproj`) contained in the repository root folder or they are run from within the **nextflow** pipelines in the root folder.

The scripts to create all figures are contained in the sub-folder `./fig` and are executable R-scripts that should be run from the root folder of the project, eg:

```sh
Rscript --vanilla R/fig/plot_XX.R input_parameter1 input_parameter2
```

Please refer to the [documentation](https://k-hench.github.io/hamlet_radiation/) for some more in depth explanations.

---

<p align="center"><img src="../logo.svg" alt="logo" width="150"/></p>