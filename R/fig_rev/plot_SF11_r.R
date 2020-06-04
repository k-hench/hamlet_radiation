#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF11.R 2_analysis/newhyb/nh_input/NH.Results/
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/newhyb/nh_input/NH.Results/")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(prismatic)
library(paletteer)
library(patchwork)

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

# p_all <- p_loc[[1]] + 
#   p_loc[[2]] + 
#   p_loc[[3]] +
#   guide_area() +
#   plot_layout(ncol = 1,
#               heights = c(3,6,1.5,.5),
#               guides = "collect")

p <- (((p_loc[[1]]  + theme(legend.position = "none")) +
    (p_loc[[3]]  + theme(legend.position = "none")) +
    plot_layout(ncol = 1,heights = c(1,.4))) |
  (p_loc[[2]] + theme(legend.position = c(.025,.0125), legend.justification = c(0,0)))) + 
  plot_annotation(tag_levels = 'a')

ggsave(filename = "figures/SF11.pdf",
       plot = p,
       height = 16,
       width = 22, 
       device = cairo_pdf)
