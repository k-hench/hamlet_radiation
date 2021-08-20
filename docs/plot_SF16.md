---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 16



## Summary

This is the accessory documentation of Figure S16.
The Figure can be recreated by running the **R** script `plot_SF16.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF16.R 2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support
```

## Details of `plot_SF16.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen), [**ape**](http://ape-package.ird.fr/) and [**ggtree**](https://github.com/YuLab-SMU/ggtree).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF16.R 2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support
# ===============================================================
# This script produces Suppl. Figure 16 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support")
# script_name <- "R/fig/plot_SF16.R"
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF16.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support
#> ────────────────────────────────────────── /current/working/directory ───
```

The path of the tree file is received and stored inside a more descriptive variable.


```r
# config -----------------------
tree_serr_file <- as.character(args[1])
```

Then, the default tree layout is defined and the default tree color is set.


```r
clr_neutral <- rgb(.6, .6, .6)
lyout <- 'circular'
```

Then, the tree file is read and the tree is rooted once with the Serranus samples as outgroup and once midpoint-rooted.


```r
tree_s <- read.tree(tree_serr_file)
tree_s_rooted <- root(tree_s, outgroup = c("28393torpan", "s_tort_3torpan", "20478tabhon" ))
tree_s_mid <- phangorn::midpoint(tree_s_rooted)
```

Now, the support values of the tree are transformed into discrete support classes.


```r
tree_s_data <- open_tree(ggtree(tree_s_mid, layout = "circular"), 180) %>% 
  .$data %>% 
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         support = as.numeric(label),
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% 
           factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")))
```

At this point, the basic pylogenetic tree can be drawn.


```r
p_s_tree <- rotate_tree(open_tree(ggtree(tree_s_data,
                                          aes(color = spec),
                                          layout = "circular"), 
                                   180), 0) +
  geom_tiplab2(aes(label = str_sub(label, -6, -1)),
               size = 3, offset = .001) +
  geom_nodepoint(aes(fill = support_class,
                     size = support_class),
                 shape = 21) +
  ggtree::geom_treescale(width = .05,
                         x = -.02,
                         y = 158, 
                         offset = -15, fontsize = 3,
                         color = clr_neutral) +
  scale_color_manual(values = c(ungrouped = clr_neutral, 
                                GenomicOriginsScripts::clr2),
                     guide = FALSE) +
  scale_fill_manual(values = c(`(0,50]` = "transparent",
                               `(50,70]` = "white",
                               `(70,90]` = "gray",
                               `(90,100]` = "black"),
                    drop = FALSE) +
  scale_size_manual(values = c(`(0,50]` = 0,
                               `(50,70]` = 1.5,
                               `(70,90]` = 1.5,
                               `(90,100]` = 1.5),
                    na.value = 0,
                    drop = FALSE) +
  guides(fill = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2),
         size = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2)) +
  theme_void()
```

Unfortunately, the polar coordinate system underlying the circular tree layout introduces a lot of empty space for a open tree that spans 180 degrees.
In order to crop that empty space, the basic tree is converted into grid object and used as an annotation in a different ggplot.
This uses cartesian coordinates and is easily cropped to create the final figure.


```r
y_sep <- .05
x_shift <- -.03
p_done <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(-.01, .54),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(p_s_tree + theme(legend.position = "none")),
                    ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
                    xmin = -.1, xmax = 1.1) +
  annotation_custom(grob = cowplot::get_legend(p_s_tree),
                    ymin = .05, ymax = .15,
                    xmin = .4, xmax = .6) +
  theme_void()
```



Finally, we can export Figure S16.


```r
scl <- 1.5
hypo_save(plot = p_done,
          filename = "figures/SF16.pdf",
          width = 7.5 * scl,
          height = 4 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
```

---
