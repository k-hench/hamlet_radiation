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

Rscript --vanilla R/fig/plot_SF1.R 2_analysis/pca/

```

## Details of `plot_SF1.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/) 

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF1.R 2_analysis/pca/
# ===============================================================
# This script produces Suppl. Figure 1 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/pca/")
# script_name <- "R/fig/plot_SF1.R"
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(patchwork)
library(hypoimg)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF1.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/pca/
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
pca_dir <- as.character(args[1])
clr_alt <- clr
clr_alt[["uni"]] <- "lightgray"
```

Then the layout of the legend elements, labels and patch sizes are pre-configured.


```r
fish_tib <- tibble(short = names(clr)[!names(clr) %in% c("flo", "tab", "tor")],
                   x = c(0.5,  3.5,  7,  9.7, 12.25, 15.25, 18, 21.5))

key_sz <- .75
sp_fam <- rep(c("H", "S", "H"), c(8, 2, 1)) %>% set_names(nm = names(sp_names))
```

Now, the legend is created.


```r
p_leg <- fish_tib %>% 
  ggplot() +
  coord_equal(xlim = c(-.05, 24), expand = 0) +
 geom_tile(aes(x = x, y = 0,
                fill = short, 
                color = after_scale(prismatic::clr_darken(fill, .25))),
            width = key_sz, height = key_sz, size = .3) +
  geom_text(aes(x = x + .6, y = 0,
                label = str_c(sp_fam[short], ". ", sp_names[short])), 
            hjust = 0, fontface = "italic", size = plot_text_size / ggplot2:::.pt) +
  pmap(fish_tib, plot_fish_lwd, width = 1, height = 1, y = 0) +
  scale_fill_manual(values = clr, guide = FALSE) +
  theme_void()
```



Then, by running the function `pca_plot_no_fish()` one on each location, the figure can be assembled.


```r
p_done <- cowplot::plot_grid((tibble(loc = c("bel.", "hon.", "pan."), 
                              mode = rep(c("subset_non_diverged"), 3),
                              pc1 = 1,
                              pc2 = 2) %>% 
                        pmap(pca_plot_no_fish) %>% 
                        wrap_plots(ncol = 3) +
                        plot_annotation(tag_levels = "a") & 
                        theme(plot.background = element_blank())),
                     p_leg,
  ncol = 1 ,
  rel_heights = c(1,.1))
```



Finally, we can export Figure S1.


```r
hypo_save(p_done, filename = 'figures/SF1.pdf',
          width = f_width,
          height = f_width * .38,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
```
