#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla network_generate_background_graphs.R ${n_snps} \$BASE_DIR/R/network_functions.R \$BASE_DIR/R/project_config.R ${vcf} sample_file.txt
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c(50,'network_test/network_functions.R','~/Desktop/chapter2/R/project_config.R', '/software/KBIN/downsamplevcf.jar', 'graph.vcf.gz', 'network_test/samples_no_out.txt')
args = commandArgs(trailingOnly=FALSE)
args = args[7:12]
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
java_downsample <- as.character(args[4])
vcf_file <- as.character(args[5])
sample_file <- as.character(args[6])
source(functions_script)
source(proj_config)
# load data -------------------

samples <- read_tsv(sample_file,col_names = 'ID') %>%
  mutate(pop = str_sub(ID,start = -3,end = -1),
         spec =  str_sub(ID,start = -6,end = -4),
         grp = str_c(spec,pop))

n_snps %>%
  purrr::map(generate_random_subset_graphs, n_graphs = 50)
