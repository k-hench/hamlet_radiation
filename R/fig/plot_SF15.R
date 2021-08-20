#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF15.R \
#     2_analysis/raxml/lg04.1_hySN.raxml.support \
#     2_analysis/raxml/lg12.3_hySN.raxml.support \
#     2_analysis/raxml/lg12.4_hySN.raxml.support
# ===============================================================
# This script produces Suppl. Figure 15 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/raxml/lg04.1_hySN.raxml.support",
#           "2_analysis/raxml/lg12.3_hySN.raxml.support",
#           "2_analysis/raxml/lg12.4_hySN.raxml.support")
# script_name <- "R/fig/plot_SF15.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
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


p1 <- plot_outl_tree_s(tree_data[[1]])
p2 <- plot_outl_tree_s(tree_data[[2]], show_legend = FALSE)
p3 <- plot_outl_tree_s(tree_data[[3]], show_legend = FALSE)

p_done <- p1 + p2 + p3 + plot_annotation(tag_levels = 'a') + plot_layout(ncol = 1)

hypo_save(plot = p_done,
          filename = "figures/SF15.pdf",
          width = f_width,
          height = f_width * 1.5,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
