#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/ccf_pairs.R
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
library(tidyverse)

targets <- read_tsv("metadata/sample_info.txt") %>%
  filter(id %in% c("18161", "18172", "18174", "18152", "18180")) %>%
  .$label

querries <- read_tsv("metadata/sample_info.txt") %>%
  filter(id %in% c("17996" ,"17998", "18195", "18226", "18237",
                   "18267" ,"18274", "18276", "20092", "20128",
                   "18426" ,"18429", "18432", "18912", "27678",
                   c("18161", "18172", "18174", "18152", "18180")))%>%
  .$label

filter <- function(x, y) x == y

cross_df(list(target = targets, querry = querries), .filter = filter) %>%
  write_tsv("ressources/plugin/ccf_pairs.tsv", col_names = TRUE)
