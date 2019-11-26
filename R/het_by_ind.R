#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla het_by_ind.R <ind>.het.gz 50000
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c("het_import/17997indbel.het.gz", '50000')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(tidyverse)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

args <- process_input(script_name, args)
# config -----------------------
het_file <- as.character(args[1])
win_sz <- as.numeric(args[2])

# load data -------------------
bin_by_position <- function(x, window_size = 1000) ((x-1)%/%window_size)+1

ind_id <- het_file %>% 
  str_remove('^.*/') %>% 
  str_remove('.het.gz')

data <- vroom::vroom(het_file, delim = '\t') %>%
  mutate(window = bin_by_position(POS, window_size = win_sz)) %>%
  group_by(CHROM,window,IND) %>%
  summarise(het = mean(HET==.5)) %>%
  ungroup() %>%
  left_join(hypogen::hypo_karyotype) %>%
  mutate(win_start = ((window-1)*win_sz)+1,
         win_end = ifelse(window*win_sz < LENGTH, window*win_sz, LENGTH),
         win_pos = (win_start + win_end)/2,
         spec = str_sub(IND, -6,-4),
         loc = str_sub(IND, -3,-1),
         win_id = str_c(CHROM,window,
                        format(win_start, scientific = FALSE, trim = TRUE),
                        format(win_end, scientific = FALSE, trim = TRUE),
                        sep = '_')) %>%
  select(win_id,het) %>%
  set_names(nm = c('win_id', ind_id))

write_tsv(x = data, path = str_c(ind_id, '_win_het.gz'))


