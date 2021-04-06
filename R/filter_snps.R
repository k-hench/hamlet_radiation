#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla filter_snps.R in.weir.fst.gz 80 out_prefix
# ===============================================================
# This script produces Figure 1 of the study "The genomic onset of a marine radiation"
# by Hench, McMillan and Puebla
#   ---------------------------------------------------------------
# ===============================================================
# args <- c('in.weir.fst.gz', "80", filterSet.80SNPs.',name,'.snps')
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

args <- process_input(script_name, args)

# config -----------------------
in_file <- as.character(args[1])
nr_snps <- as.numeric(args[2])
out_prefix <- as.character(args[3])

# start script -------------------

vroom::vroom(in_file, delim = '\t') %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>%
  mutate(FST_RANK = rank(-WEIR_AND_COCKERHAM_FST, ties.method = "random")) %>%
  select(CHROM, POS, WEIR_AND_COCKERHAM_FST, FST_RANK) %>% 
  filter(FST_RANK < nr_snps)  %>% 
  select(CHROM, POS) %>%
  write_tsv(path = str_c(out_prefix, "_", nr_snps, "SNPs",'.snps' ))