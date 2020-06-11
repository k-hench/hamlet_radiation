#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF1.R 2_analysis/newhyb/nh_input/NH.Results/
# ===============================================================
# This script produces Suppl. Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/newhyb/nh_input/NH.Results/")
# script_name <- "plot_SF1.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(prismatic)
library(paletteer)
library(patchwork)
library(ggtext)
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

# config -----------------------
base_dir <- as.character(args[1])

folders <- dir(base_dir)

p_loc <- c("bel", "hon", "pan") %>%
  map(plot_loc)

p <- (p_loc[[1]] +  guides(fill = guide_legend(title = "Hybrid Class")) + theme_hyb(legend.position = c(1,1)) ) +
  (p_loc[[2]] + theme_hyb() ) +
  (p_loc[[3]] + theme_hyb() )  +
  plot_layout(ncol = 1, heights = c(10,15,3) %>% label_spacer())+
  plot_annotation(tag_levels = 'a')

hypo_save(filename = "figures/SF1.pdf",
       plot = p,
       height = 16,
       width = 10,
       device = cairo_pdf,
       comment = plot_comment)
