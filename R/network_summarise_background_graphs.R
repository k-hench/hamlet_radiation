#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla network_summarise_background_graphs.R \$BASE_DIR/R/network_functions.R \$BASE_DIR/R/project_config.R
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c('network_test/network_functions.R','~/Desktop/chapter2/R/project_config.R')
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]
print(args)
# setup -----------------------
library(vcfR)
library(adegenet)
library(gstudio)
library(popgraph)
library(tidyverse)
library(igraph)

# config -----------------------
functions_script <- as.character(args[1])
proj_config <- as.character(args[2])
graph_dir <- as.character(args[3])
source(functions_script)
source(proj_config)
# load data -------------------
files <- dir( graph_dir, pattern = 'popgr.background\\.' )

summary_data <- files %>%
  str_c(graph_dir,.) %>%
  purrr::map(summarise_popgr) %>%
  bind_rows()

plot_data <- summary_data %>%
  gather(key = 'key',
         value = 'value',7:15) %>%
  mutate(graph_snps = str_pad(graph_snps,pad  = '0',width = 4))

plot_data_count <- plot_data %>%
  group_by(graph_snps, connected, key) %>%
  count()

p1 <- plot_data %>%
  ggplot(aes(x = connected, y= value))+
  geom_boxplot(aes(fill = factor(connected))) +
  geom_text(data = plot_data_count, aes(y = -Inf, label = n),vjust = .01 )+
  facet_grid(key ~ graph_snps, scales = 'free') +
  scale_fill_brewer(palette = 'PuOr', guide = FALSE)

#ggsave(plot = p1, filename = 'network_background.png', width = 15, height = 10, dpi = 250)
ggsave(plot = p1, filename = 'network_background.pdf', width = 15, height = 10)
write_tsv(x = summary_data, path = "network_background_data.tsv")