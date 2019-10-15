#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla table_fst_outliers.R multi_fst.50k.tsv.gz
# ===============================================================
# This script produces Fst outlier table of the study 
# "The genomic origins of a marine radiation"
# by Hench, McMillan an Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('figures/data/fst/multi_fst.50k.tsv.gz')
args = commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(tidyverse)
library(hypogen)
library(vroom)

cat('\n')

args <- process_input(script_name, args)

# config -----------------------
fst_file <- as.character(args[1])

outliers <-  vroom::vroom(fst_file,delim = '\t') %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST) %>%
  setNames(., nm = c('chrom', 'start', 'end', 'n_snps', 'fst') ) %>%
  left_join(.,hypogen::hypo_chrom_start, by = c('chrom' = 'CHROM')) %>%
  mutate(gpos = (start+end)/2 + GSTART) %>%
  group_by(chrom) %>%
  mutate(norm_fst = fst-mean(fst)) %>%
  ungroup() %>% 
  filter(norm_fst >= quantile(norm_fst, .998)) %>%
  group_by(chrom) %>%
  mutate(check = gpos > lag(gpos,default = 0) + 50000,
         id = cumsum(check),
         gid = str_c(chrom,'_',id)) %>%
  group_by(gid) %>%
  summarise(chrom = chrom[1],
            start = min(start),
            end = max(end),
            gstart = min(start)+GSTART[1],
            gend = max(end)+GSTART[1]) %>%
  mutate(gpos = (gstart+gend)/2)


write_tsv(outliers, 'fst_outliers_998.tsv')
