---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 11) Visualization

After all `nextflow` pipelines are successfully run to completion, each Figure (and Suppl. Figure) of the manuscript can be re-created with its respective `R` script located under `R/fig`.
These are executable `R` scripts that can be launched from the base directory;

```sh
Rscript --vanilla R/fig/plot_Fxyz.R input1 input2 ...
```

For convenience, there also exists a `bash` script that can be used to re-create all Figures in one go (git 11):

```sh
cd $BASE_DIR
bash sh/create_figures.sh
```

After running `create_figures.sh`, Figures 1 - 5 and Suppl. Figures 1 - 9 should be created withing the folder `figures/`.

In the remaining documentation, the individual Visualization scripts are going to discussed in detail.

---
