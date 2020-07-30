#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla major_allele.R aa/allele_counts.tsv aa/phased3_annotations.bed
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c("aa/allele_counts.tsv","aa/phased3_annotations.bed")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
args <- args[7:length(args)]

library(tidyverse)
# config -----------------------
allele_counts <- as.character(args[1])
bed <- as.character(args[2])

data_allele <- read_tsv(allele_counts) %>%
  separate(col = ALLELE1, into =  c("A1", "FREQ1"), convert = TRUE, sep = ":" ) %>%
  separate(col = ALLELE2, into =  c("A2", "FREQ2"), convert = TRUE, sep = ":" ) %>%
  mutate(majA = ifelse(FREQ1 > FREQ2, A1,A2)) %>%
  select(POS,majA) %>%
  set_names(nm = c("TO", "majA"))

read_tsv(bed) %>%
  left_join(.,data_allele) %>%
  mutate(AA = ifelse(AA == ".", majA, AA)) %>%
  select(-majA) %>%
  write_tsv(path = str_replace(bed,".bed$","_maj.bed"), 
            col_names = TRUE)
