---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 9



## Summary

This is the accessory documentation of Supplementary Figure .
The Figure can be recreated by running the **R** script `plot_SF.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF

```

## Details of `plot_SF.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:



The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) is loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.



```r
#> ── Script: scripts/plot_SF.R ────────────────────────────────────────────
#> Parameters read:
#>
#> ─────────────────────────────────────────── /current/working/directory ──
```




---
