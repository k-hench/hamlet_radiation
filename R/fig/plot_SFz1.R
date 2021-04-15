#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFz1.R 2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support
# ===============================================================
# This script produces Suppl. Figure 9 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/raxml/hyS_n_0.33_mac4_5kb.raxml.support")
# script_name <- "R/fig/plot_SF9.R"
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

# config -----------------------
tree_serr_file <- as.character(args[1])

clr_neutral <- rgb(.6, .6, .6)
lyout <- 'circular'

lab2spec <- function(label){
  x <- str_sub(label, start = -6, end = -4) %>% str_remove(.,"[0-9.]{1,3}$") %>% str_remove(.," ")
  ifelse(x == "",'ungrouped', x)
}

tree_s <- read.tree(tree_serr_file)
tree_s_rooted <- root(tree_s, outgroup = c("28393torpan", "s_tort_3torpan", "20478tabhon" ))
# tree_s_mid <- phangorn::midtree_s_rootedpoint(tree_s_rooted)

tree_s_data <- open_tree(ggtree(tree_s_rooted, layout = "circular"), 180) %>% 
  .$data %>% 
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         support = as.numeric(label),
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% 
           factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")))

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

y_sep <- .05
x_shift <- -.03
p_single <- ggplot() +
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

scl <- 1.5
hypo_save(plot = p_single,
          filename = "figures/SFz1.pdf",
          width = 7.5 * scl,
          height = 4 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)