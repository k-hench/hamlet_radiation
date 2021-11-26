---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 13






## Summary

This is the accessory documentation of Figure S13.
The Figure can be recreated by running the **R** script `plot_SF13.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF13.R \
    2_analysis/summaries/fst_outliers_998.tsv
```

## Details of `plot_SF13.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF13.R \
#     2_analysis/summaries/fst_outliers_998.tsv
# ===============================================================
# This script produces Suppl. Figure 13 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/summaries/fst_outliers_998.tsv")
# script_name <- "R/fig/plot_SF13.R"
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
library(tidygraph)
library(ggraph)
library(prismatic)
library(patchwork)
library(IRanges)
library(plyranges)
library(hypogen)
library(hypoimg)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF13.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/summaries/fst_outliers_998.tsv
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
outlier_file <- as.character(args[1])
outlier_regions <- read_tsv(outlier_file)
```



```r
hap_to_perc <- 100 / (166 * 165)
iterations <- c("25/10 kb","10/5 kb","15/7.5 kb") %>% set_names(value = c("7", "8", "10"))
idx <- 10
```



```r
import_map1 <- function(idx, filtmode = "bed95"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>%
    mutate(ibd_total = (ibd2_cM_m1 + 0.5*ibd1_cM_m1) / (ibd0_cM_m1 + ibd1_cM_m1 + ibd2_cM_m1))
}

import_map2 <- function(idx, filtmode = "bed95"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>%
    mutate(ibd_total = (ibd2_cM_m2 + 0.5*ibd1_cM_m2) / (ibd0_cM_m2 + ibd1_cM_m2 + ibd2_cM_m2))
}

import_bp <- function(idx, filtmode = "bed95"){
    read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_summary.tsv")) %>%
    mutate(ibd_total = (ibd2_bp + 0.5*ibd1_bp) / (ibd0_bp + ibd1_bp + ibd2_bp))
}

import_truffle <- function(idx, filtmode = "direct"){
  itteration_names <- c(str_c("10-",6:3),"7","8","9","10")
  read_tsv(glue::glue("2_analysis/ibd/no_outgr_{filtmode}_{itteration_names[idx]}.ibd.tsv")) %>%
    mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2))
}
```



```r
iterations <- c(str_c("2/5*10^",6:3," BP"), "25/10 kb", "10/5 kb", "7-5/3 kb", "15/7.5 kb")

plot_network <- function(idx, filt = 0, import_fun = import_map1,
                         x = "cM_map1", filtmode = "direct",
                         x_ax = TRUE, y_ax = TRUE, ...){
  clr2 <- GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))]
  clr2["uni"] <- rgb(.9,.9,.9)

  data <- import_fun(idx, filtmode = filtmode)

  set.seed(42)

  p <- data %>%
    as_tbl_graph() %E>%
    filter(ibd_total > filt) %N>%
    mutate(spec = str_sub(name,-6,-4),
           loc = str_sub(name,-3,-1))  %>%
    ggraph( layout = 'fr', weights = ibd_total) +
    geom_edge_link(aes(alpha = ibd_total), color = rgb(.1,.1,.1), edge_width = .15) +
    geom_node_point(aes(fill = spec,
                        shape = loc, color = after_scale(clr_darken(fill,.3))), size = .7) +
    labs(y = glue::glue("Seq. Length: {iterations[idx]}"),
         x = x) +
    scale_fill_manual("Species", values = GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))],
                      labels = GenomicOriginsScripts::sp_labs)+
    scale_edge_alpha_continuous(range = c(0,1), guide = "none") +
    scale_shape_manual("Site", values = 21:23, labels = GenomicOriginsScripts::loc_names) +
    scale_x_continuous(position = "top") +
    guides(fill = guide_legend(title.position = "top",
                               nrow = 2, override.aes = list(shape = 21, size = 2.5)),
           shape = guide_legend(title.position = "top",
                                nrow = 2, override.aes = list(size = 2.5))) +
    coord_equal()  +
    theme(panel.background = element_blank(),
          axis.title.y = element_text(),
          axis.title.x = element_text())

  if(!x_ax){ p <- p + theme(axis.title.x = element_blank())}
  if(!y_ax){ p <- p + theme(axis.title.y = element_blank())}
  p
}
```



```r
plts_cM <- tibble(idx = rep(c(7, 10, 8), each = 3),
       import_fun = rep(list(import_bp, import_map1, import_map2), 3),
       x = rep(c("bp_cM_filt.", "cM_map1", "cM_map2"), 3),
       x_ax = rep(c(TRUE, FALSE), c(3, 6)),
       y_ax = rep(FALSE, 9),
       filtmode = "bed95") %>%
  bind_rows(tibble(idx = rep(c(5, 8, 6), 2),
                   import_fun = rep(list(import_truffle), 6),
                   x = rep(c("truffle", "bed95"), each = 3),
                   x_ax = rep(rep(c(TRUE, FALSE), 1:2), 2),
                   y_ax = rep(c(TRUE, FALSE), each = 3),
                   filtmode = rep(c("direct", "bed95"), each = 3)) ) %>%
  left_join(tibble(x = c("truffle", "bed95", "bp_cM_filt.", "cM_map1", "cM_map2"),
            plot_order = seq_along(x))) %>%
  arrange(plot_order) %>%
  pmap(plot_network)
```



```r
p_done <- plts_cM %>%
  wrap_plots(nrow = 3,
             byrow = FALSE,
             guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(size = plot_text_size),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.text.align = 0,
        plot.subtitle = element_text())
```

Finally, we can export Figure S13.


```r
hypo_save(plot = p_done,
          filename = "figures/SF13.png",
          width = f_width,
          height = .75*f_width,
          dpi = 600,
          type = "cairo",
          bg = "transparent",
          comment = plot_comment)

system("convert figures/SF13.png figures/SF13.pdf")
system("rm figures/SF13.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SF13.pdf")
system(create_metadata)
```
