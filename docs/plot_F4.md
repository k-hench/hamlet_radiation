---
output: html_document
editor_options:
  chunk_output_type: console
---
# Figure 4






## Summary

This is the accessory documentation of Figure 4.

The Figure can be recreated by running the **R** script `plot_F4.R` from a (`bash` terminal):

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_F4.R \
    2_analysis/astral/astral_5000x_5kb_v1_noS.tre \
    2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv
```

## Details of `plot_F4.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), [**BAMMtools**](https://cran.r-project.org/web/packages/BAMMtools/) and on the package [**hypoimg**](https://k-hench.github.io/hypoimg).

### Config

The scripts start with a header that contains copy & paste templates to execute interactively or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F4.R \
#     2_analysis/astral/astral_5000x_5kb_v1_noS.tre \
#     2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv
# ===============================================================
# This script produces Figure 4 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/astral/astral_5000x_5kb_v1_noS.tre",
#           "2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv")
# script_name <- "R/fig/plot_F4.R"
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
#> ── Script: R/fig/plot_F4.R ──────────────────────────────────────────────
#> Parameters read:
#>  ★ 1: 2_analysis/astral/astral_5000x_5kb_v1_noS.tre
#>  ★ 2: 2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv
#> ─────────────────────────────────────────── /current/working/directory ──
```
The directories for the different data types are received and stored in respective variables.
Also, we set a few parameters for the plot layout:


```r
# config -----------------------
tree_hypo_file <- as.character(args[1])
ibd_file <- as.character(args[2])
```

### Actual Script Start


```r
tree <- read.tree(tree_hypo_file)
tree$edge.length <- replace(tree$edge.length, tree$edge.length == "NaN", 0.05)   # Set terminal branches to 0.05
```



```r
tree_rooted <- root(phy = tree, outgroup = "PL17_160floflo")
clr_neutral <- rgb(.2, .2, .2)
```



```r
### Prepare tree and categorize support values
tree_plus <- ggtree(tree_rooted) %>%
  ggtree::rotate(node = 205) %>%
  .$data %>%
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
         support = as.numeric(label) * 100,
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
  )
```



```r
t_plot <- (ggtree(tree_plus, layout = 'fan', aes(color = spec), size = .2) %>% 
             ggtree::rotate_tree(angle = -100)) +
  geom_tippoint(aes(color = spec,
                    shape = loc,
                    fill = after_scale(color)), size = .5) +
  geom_nodepoint(data = tree_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
                 aes(fill = support_class,
                     size = support_class),
                 shape = 21,
                 color = clr_neutral) +
  scale_color_manual(values = c(GenomicOriginsScripts::clr2, ungrouped = "gray60"),
                     guide = 'none') +
  scale_shape_manual(values = c(bel = 21, flo = 24, hon = 22, pan = 23), labels = GenomicOriginsScripts::loc_names,
                     guide = 'none') +
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
  # Add scale bar:
  ggtree::geom_treescale(width = .05,
                         x = .13, y = 130,
                         offset = -1,
                         linesize = .2,
                         fontsize = plot_text_size/.pt,
                         color = clr_neutral) +
  scale_x_continuous(limits = c(-.05, .26), expand = c(0,0)) +
  guides(fill = guide_legend(title = "Node Support Class", title.position = "top",
                             nrow = 2),
         size = guide_legend(title = "Node Support Class", title.position = "top",
                             nrow = 2)) +
  theme_void() +
  theme(legend.title.align = 0.5,
        legend.text = element_text(color = "gray20"),
        legend.title = element_text(color = "gray20"))
```



```r
y_sep <- .55
x_shift <- .5
```



```r
p1 <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
                    ymin = -.15 - y_sep , 
                    ymax = 1 + y_sep,
                    xmin = .25 - x_shift,
                    xmax = 1 + x_shift) +
  theme_void()
```



```r
data_ibd <- read_tsv(ibd_file) %>% 
  mutate(ibd_total = (ibd2_cM_m1 + 0.5*ibd1_cM_m1) / (ibd0_cM_m1 + ibd1_cM_m1 + ibd2_cM_m1))
```



```r
set.seed(42)
p2 <- data_ibd %>% 
  as_tbl_graph() %>%
  mutate(spec = factor(str_sub(name, -6, -4), levels = names(clr)[-10:-9]), 
         loc = factor(str_sub(name, -3, -1), levels = c("bel", "hon", "pan", "flo")))  %>% 
  ggraph( layout = 'fr', weights = ibd_total) +
  geom_edge_link(aes(alpha = ibd_total, edge_width = ibd_total), color = rgb(.1,.1,.1)) +
  geom_node_point(aes(fill = spec,
                      shape = loc, color = after_scale(clr_darken(fill,.3))), size = 1.2) +
  scale_fill_manual("Species", values = GenomicOriginsScripts::clr2[!(names(GenomicOriginsScripts::clr2) %in% c( "tor", "tab"))],
                    labels = GenomicOriginsScripts::sp_labs, drop = FALSE)+
  scale_edge_alpha_continuous(range = c(0, 1), guide = "none") +
  scale_edge_width_continuous(range = c(.1, .4), guide = "none") +
  scale_shape_manual("Site", values = c(bel = 21, flo = 24, hon = 22, pan = 23),
                     labels = GenomicOriginsScripts::loc_names, drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 21, size = 2.5), title.position = "top",label.hjust = 0),
         shape = guide_legend(nrow = 2, title.position = "top")) +
  coord_fixed() +
  theme(text = element_text(size = GenomicOriginsScripts::plot_text_size),
        panel.background = element_blank())
```



```r
p_done <- cowplot::plot_grid(cowplot::plot_grid(p1, p2 + theme(legend.position = "none"), rel_widths = c(1,.9),
                                                labels = c("a", "b"), label_fontface = "plain", label_size = plot_text_size),
                             cowplot::plot_grid(cowplot::get_legend(t_plot +
                                                                      theme_minimal(base_size = GenomicOriginsScripts::plot_text_size) +
                                                                      theme(legend.key.width = unit(3,"pt"))),
                                                cowplot::get_legend(p2 +
                                                                      theme_minimal(base_size = GenomicOriginsScripts::plot_text_size) +
                                                                      theme(legend.position = "bottom",
                                                                            legend.key.width = unit(3,"pt"))),
                                                rel_widths = c(.3,1)),
                             rel_heights = c(1,.15),
                             ncol = 1)
```

Finally, we can export Figure 4.


```r
hypo_save(plot = p_done,
          filename = "figures/F4.png",
          width = GenomicOriginsScripts::f_width,
          height = GenomicOriginsScripts::f_width * .65,
          bg = "white",
          type = "cairo",
          dpi = 600,
          comment = plot_comment)

system("convert figures/F4.png figures/F4.pdf")
system("rm figures/F4.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/F4.pdf")
system(create_metadata)
```

---
