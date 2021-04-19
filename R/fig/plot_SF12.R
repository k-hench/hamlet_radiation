#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF12.R 2_analysis/raxml/lg04.1_hySN.raxml.support \
#    2_analysis/raxml/lg12.3_hySN.raxml.support \
#    2_analysis/raxml/lg12.4_hySN.raxml.support
# ===============================================================
# This script produces Suppl. Figure 12 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/raxml/lg04.1_hySN.raxml.support",
#           "2_analysis/raxml/lg12.3_hySN.raxml.support",
#           "2_analysis/raxml/lg12.4_hySN.raxml.support")
# script_name <- "R/fig/plot_SF12.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)
library(patchwork)
library(phangorn)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
tree_file_lg04_1 <- as.character(args[1])
tree_file_lg12_3 <- as.character(args[2])
tree_file_lg12_4 <- as.character(args[3])

trees <- c(tree_file_lg04_1, tree_file_lg12_3, tree_file_lg12_4) %>% 
  map(.f = function(file){
    read.tree(file) %>%
      root(phy = ., outgroup = c("28393torpan", "s_tort_3torpan", "20478tabhon" )) %>% 
    midpoint()}
  )

clr_neutral <- rgb(.6, .6, .6)
lyout <- 'circular'

lab2spec <- function(label){
  x <- str_sub(label, start = -6, end = -4) %>% str_remove(.,"[0-9.]{1,3}$") %>% str_remove(.," ")
  ifelse(x == "",'ungrouped', x)
}

tree_data <- trees %>% 
  map(.f = function(tree_in){
    open_tree(ggtree(tree_in, 
                     layout = lyout), 180) %>%
      .$data %>% 
      mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
             support = as.numeric(label),
             support_class = cut(support, c(0,50,70,90,100)) %>% 
               as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
      )}
  )

plot_outl_tree <- function(tree, show_legend = TRUE){
  p <- open_tree(ggtree(tree,
              layout = lyout,
              aes(color = spec), size = .3), 180) +
    geom_tiplab(aes(color = lab2spec(label), 
                    label = str_sub(label, -6, -1)),
                size = 1.7,
                hjust = -.1)+
    ggtree::geom_treescale(width = .01,
                           x = .001, 
                           y = 158, 
                           offset = -1,
                           fontsize = 1.5,
                           color = clr_neutral) +
    ggtree::geom_nodepoint(aes(fill = support_class, size = support_class),
                           shape = 21) +
    scale_color_manual(values = c(ungrouped = clr_neutral, 
                                  GenomicOriginsScripts::clr2),
                       guide = FALSE) +
    scale_fill_manual(values = c(`(0,50]` = "transparent",
                                 `(50,70]` = "white",
                                 `(70,90]` = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE) +
    scale_size_manual(values = c(`(0,50]` = 0,
                                 `(50,70]` = 1,
                                 `(70,90]` = 1,
                                 `(90,100]` = 1),
                      na.value = 0,
                      drop = FALSE) +
    guides(fill = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2),
           size = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2)) +
    theme_void() 
    # theme_bw() +
  
  if(show_legend){
    p <- p +
      theme(legend.position = c(.5, .75),
            legend.justification = c(.5, 1),
            legend.text = element_text(size = plot_text_size))
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  ggplot() +
    coord_equal(xlim = c(0, 1),
                ylim = c(-.01, .54),
                expand = 0) +
    annotation_custom(grob = ggplotGrob(p),
                      ymin = -.6, ymax = .6,
                      xmin = -.1, xmax = 1.1) +
    theme_void()
}

p1 <- plot_outl_tree(tree_data[[1]])
p2 <- plot_outl_tree(tree_data[[2]], show_legend = FALSE)
p3 <- plot_outl_tree(tree_data[[3]], show_legend = FALSE)

p_done <- p1 + p2 + p3 + plot_annotation(tag_levels = 'a') + plot_layout(ncol = 1)

scl <- 2
hypo_save(plot = p_done,
          filename = "figures/SF12.pdf",
          width = f_width_half * scl,
          height = f_width_half * 1.5 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)