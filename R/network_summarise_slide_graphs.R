#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanillanetwork_summarise_slide_graphs.R ${n_snps} \$BASE_DIR/R/network_functions.R \$BASE_DIR/R/project_config.R window_file.nr_sps.tsv
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c(50,'network_test/network_functions.R','~/Desktop/chapter2/R/project_config.R', 'network_test/LG01_quarter.snp_windows.500.tsv.gz')
args = commandArgs(trailingOnly=FALSE)
args = args[7:10]
print(args)
# setup -----------------------
library(vcfR)
library(adegenet)
library(gstudio)
library(popgraph)
library(tidyverse)
library(hypogen)
library(igraph)

# config -----------------------
n_snps <- as.numeric(args[1])
functions_script <- as.character(args[2])
proj_config <- as.character(args[3])
window_file <- as.character(args[4])

source(functions_script)
source(proj_config)
# -----------------------------
network_windows <-  read_tsv(window_file)

windows_pos <- network_windows %>%
  select(lg, gwin_id, start, end) %>%
  left_join(hypogen::hypo_chrom_start,by = c('lg' = 'CHROM')) %>%
  mutate(gwin_id = as.character(gwin_id),
         POS = (start + end)/2,
         GPOS = POS + GSTART)

slide_files <- dir(pattern = 'popgr.LG')

summary_data <- slide_files %>%
  purrr::map(summarise_slide_gr) %>%
  bind_rows() %>%
  left_join(windows_pos) %>%
  gather('key','value',nr_cluster:eigenCent)

p1 <- summary_data %>%
  ggplot(aes(x=GPOS, y = value))+
  geom_hypo_LG()+
  geom_line() +
  ggtitle(str_c(n_snps,' SNPs'))+
  geom_point(aes(fill = connected),shape = 21)+
  facet_grid(key ~ ., scales = 'free') +
  scale_x_hypo_LG()+
  scale_fill_hypo_LG_bg() +
  scale_fill_brewer(palette = 'PuOr') +
  theme_hypo()

ggsave(plot = p1, filename = str_c("network_slide.", n_snps, ".png"), width = 15, height = 10, dpi = 250)
write_tsv(x = summary_data, path = str_c("network_slide_data.", n_snps, ".tsv"))