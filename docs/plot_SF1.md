---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 1






## Summary

This is the accessory documentation of Figure S1.
The Figure can be recreated by running the **R** script `plot_SF1.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF1.R \
    ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv
```

## Details of `plot_SF1.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF1.R \
#     ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv
# ===============================================================
# This script produces Suppl. Figure 1 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv" )
# script_name <- "R/fig/plot_SF1.R"
args <- commandArgs(trailingOnly = FALSE)
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name)))
args = args[7:length(args)]
cat(' ')
cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args),': ',crayon::green(args),'\n'))
cli::rule(right = getwd())
```

```r
#> ── Script: R/fig/plot_SF1.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
rate_file <- as.character(args[1])
```



```r
### Tip-specific speciation rates across Fish Tree of Life
### ------------------------------------------------------
## Import FToL data
rates <- read_csv(file = rate_file)

p_done <- ggplot(rates, aes(lambda.tv)) +
  geom_histogram(binwidth = 0.25, colour="grey20", fill="grey80", size = .2) +
  labs(x = "Mean speciation rate", y = "Number of species") +
  geom_segment(aes(x = 2.5 , y = 200, xend = 2.5, yend = 1500),
               size = 0.2, color = "#0976BA") +
  annotate("text", x = 2.5, y = 1750, label = "Hamlets",
           size =  plot_text_size / ggplot2:::.pt, color = "#0976BA") +
  geom_segment(aes(x = 3.5 , y = 200, xend = 3.5, yend = 1500),
               size = 0.2, color = "grey60") +
  annotate("text", x = 3.5, y = 2150, label = "Haplo-\nchromines",
           size =  plot_text_size / ggplot2:::.pt, color = "grey60") +
  geom_segment(aes(x = 4.5 , y = 200, xend = 4.5, yend = 1500),
               size = 0.2, color = "grey60") +
  annotate("text", x = 4.5, y = 1750, label = "Labeobarbus",
           size = plot_text_size / ggplot2:::.pt, color = "grey60") +
  coord_cartesian(xlim = c(-.2, 5.1),
                  ylim = c(-200, 8400),
                  expand = 0) +
  theme_bw(base_size = plot_text_size) +
  theme(panel.grid.minor = element_blank())
```

Finally, we can export Figure S1.


```r
hypo_save("figures/SF1.pdf",
       plot = p_done,
       width = f_width_half,
       height = f_width_half * .6,
       comment = plot_comment,
       device = cairo_pdf,
       bg = "transparent")
```
