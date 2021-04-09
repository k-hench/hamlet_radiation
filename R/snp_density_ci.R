#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/snp_density_ci.R window_stats.tsv.gz random_1k_windows.bed.gz
# ===============================================================
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/revPoMo/window_stats.tsv.gz", "random_1k_windows.bed.gz")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(tidyverse)
library(hypogen)

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
window_stats_file <- as.character(args[1])
out_name <- as.character(args[2])

data <- vroom::vroom(window_stats_file, delim = "\t") %>%
  left_join(hypo_chrom_start) %>%
  mutate(GPOS = GSTART + (START + END) / 2 )

cov_tres <- quantile(data$REL_COV, probs = .66)
snp_tres <- quantile(data$SNP_density, probs = .66)

set.seed(42)
random_subset <- data %>%
  sample_n(1000) 

random_subset %>%
  select(CHROM:END) %>%
  write_tsv(file = out_name)
