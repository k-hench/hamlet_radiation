#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/tab/table_s4.R 2_analysis/summaries/fst_outliers_998.tsv
# ===============================================================
# This script produces Table X1 of the study "The genomic origins of a marine radiation"
# by Hench, McMillan an Puebla
#   ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/summaries/fst_outliers_998.tsv')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)

# config -----------------------
fst_outliers <- as.character(args[1])
# start script -------------------

outlier_data <- read_tsv(fst_outliers)

table_out <- outlier_data %>%
  pmap_dfr(get_genes) %>% 
  ungroup() %>%
  filter(!duplicated(label)) %>%
  group_by(gid) %>%
  summarise(n_genes = length(label),
            genes = str_c(label,collapse = ', ') %>%
              str_replace_all(pattern = '([a-z\\.0-9\\{\\}\\\\]*, [a-z\\.0-9\\{\\}\\\\]*, [a-z\\.0-9\\{\\}\\\\]*, [a-z\\.0-9\\{\\}\\\\]*), ',replacement = '\\1@')) %>%
  ungroup() %>%
  left_join(outlier_data) %>%
  separate(genes, sep = '@', paste0("genes_", seq_len(8)), fill = "right") %>%
  pivot_longer(cols = genes_1:genes_8,values_to = 'genes') %>%
  filter(!(is.na(genes))) %>%
  mutate(gid = str_replace(string = gid, pattern = '_',replacement = '\\\\_')) %>%
  group_by(gid) %>%
  mutate(n_gid = length(gid),
         check = duplicated(gid)) %>%
  ungroup() %>%
  mutate_at(.vars = vars(gid, chrom, start, end, n_genes),
            .funs = list(~ ifelse(check, '', as_multirow(., n_gid)))) %>%
  select(chrom, gid, start, end, n_genes, genes) %>%
  rename(LG = chrom, Id = gid, n = n_genes, Start = start, End = end, Genes = genes)

export_2_latex(table = table_out, name = 'tables/suppl_tab4.tex')
