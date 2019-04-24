#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla network_generate_slide_graphs.R ${n_snps} \$BASE_DIR/R/network_functions.R \$BASE_DIR/R/project_config.R ${lg} ${start} ${end} ${win_id} ${gwin_id}
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c(50,'network_test/network_functions.R','~/Desktop/chapter2/R/project_config.R')
args = commandArgs(trailingOnly=FALSE)
args = args[7:15]
print(args)
# setup -----------------------
library(vcfR)
library(adegenet)
library(gstudio)
library(popgraph)
library(tidyverse)
library(igraph)

# config -----------------------
n_snps <- as.numeric(args[1])
functions_script <- as.character(args[2])
proj_config <- as.character(args[3])
lg <- as.character(args[4])
start <- as.numeric(args[5])
end <- as.numeric(args[6])
win_id <- as.numeric(args[7])
gwin_id <- as.numeric(args[8])
sample_file <- as.character(args[9])

source(functions_script)
source(proj_config)
# -----------------------------
samples <- read_tsv(sample_file,col_names = 'ID') %>%
  mutate(pop = str_sub(ID,start = -3,end = -1),
         spec =  str_sub(ID,start = -6,end = -4),
         grp = str_c(spec,pop))

network_windows <-  tibble( lg = lg,
                            start = start,
                            end = end,
                            n_snps = n_snps,
                            win_id = win_id,
                            gwin_id = gwin_id )

network_windows %>%
  purrr::pmap(slide_graph)