---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 17






## Summary

This is the accessory documentation of Figure S17.
The Figure can be recreated by running the **R** script `plot_SF17.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF17.R \
    2_analysis/astral/astral_5000x_5kb_v1_all.tre
```

## Details of `plot_SF17.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF17.R \
#     2_analysis/astral/astral_5000x_5kb_v1_all.tre
# ===============================================================
# This script produces Suppl. Figure 17 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/astral/astral_5000x_5kb_v1_all.tre")
# script_name <- "R/fig/plot_SF17.R"
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
library(hypogen)
library(ape)
library(ggtree)
library(tidygraph)
library(ggraph)
library(patchwork)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF17.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/astral/astral_5000x_5kb_v1_all.tre
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
tree_hypo_file <- as.character(args[1])
tree <- read.tree(tree_hypo_file)
tree$edge.length <- replace(tree$edge.length, tree$edge.length == "NaN", 0.05)   # Set terminal branches to 0.05
tree$edge.length[c( 81, 83)] <- tree$edge.length[c( 81, 83)] * 0.1
```



```r
tree_rooted <- root(phy = tree, outgroup = c("s_tort_3torpan", "20478tabhon", "28393torpan"))
clr_neutral <- rgb(.2, .2, .2)
```



```r
### Prepare tree and categorize support values
tree_plus <- ggtree(tree_rooted)  %>%
  .$data %>%
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
         support = as.numeric(label) * 100,
         support_class = cut(support, c(0,50,70,90,100)) %>%
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
          `branch.length` = if_else(node %in% c( 212, 213), `branch.length` * .0001, `branch.length`),
          branch_type = if_else(node %in% c(212, 213), "broken", "whole")
         )
```



```r
t_plot <- (ggtree(tr = tree_plus,
                   layout = 'fan',
                      aes(color = spec, linetype = branch_type), size = .2)) + #%>%
    geom_tippoint(aes(color = spec,
                      shape = loc,
                      fill = after_scale(color)), size = .5) +
    geom_nodepoint(data = tree_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = clr_neutral) +
    scale_color_manual(values = c(GenomicOriginsScripts::clr2, ungrouped = "gray60"), labels = GenomicOriginsScripts::sp_labs) +
    scale_shape_manual(values = c(bel = 21, flo = 24, hon = 22, pan = 23), labels = GenomicOriginsScripts::loc_names) +
    scale_fill_manual(values = c(`(0,50]`   = "transparent",
                                 `(50,70]`  = "white",
                                 `(70,90]`  = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE) +
    scale_size_manual(values = c(`(0,50]`   = 0,
                                 `(50,70]`  = .8,
                                 `(70,90]`  = .8,
                                 `(90,100]` = .8),
                      na.value = 0,
                      drop = FALSE) +
    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
    # Add scale bar:
    ggtree::geom_treescale(width = .2,
                           x = .13, y = 85.5,
                           offset = -7,
                           linesize = .2,
                           fontsize = plot_text_size/.pt,
                           color = clr_neutral) +
    guides(fill = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2, label.hjust = 0),
           size = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2, label.hjust = 0),
           shape = guide_legend(title = "Location", title.position = "top",
                               nrow = 2, label.hjust = 0),
           color = guide_legend(title = "Species", title.position = "top",
                               ncol = 2, label.hjust = 0)) +
    theme_void() +
    theme(legend.position = 'bottom',
      legend.title.align = 0,
      legend.text = element_text(color = "gray20"),
      legend.title = element_text(color = "gray20"))
```



```r
y_sep <- .1
x_shift <- .1
p_tdone <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
                    ymin = 0 - y_sep ,
                    ymax = 1 + y_sep,
                    xmin = 0 - x_shift,
                    xmax = 1 + x_shift) +
  theme_void()
```



```r
p_done <- cowplot::plot_grid(p_tdone, cowplot::get_legend(t_plot +
                                                  theme_minimal(base_size = plot_text_size) +
                                                  theme(legend.position = "right",
                                                        legend.title.align =0,
                                                        legend.key.height = unit(7,"pt"))),rel_widths = c(1,.5))
```

Finally, we can export Figure S17.


```r
scl <- .75
hypo_save(p_done, filename = 'figures/SF17.pdf',
          width = f_width,
          height = .6 * f_width,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
```
