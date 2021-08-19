#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla dstats.R

### Extract significant Dmin, BBAA trios
### -------------------------------------

### Jul 12, 2021, Martin Helmkampf

library(tidyverse)
library(gtools)
dataset <- "ld05"
# Alternatives: 5kb, ld06, ld05

correction <- "bonferroni"
# Alternatives: bonferroni, fdr

## D-stats for topologies minimizing D (Dmin)
dmin <- read_tsv(paste("hyp_", dataset, "_dtrios_Dmin.txt", sep ="")) %>%
  rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
  mutate(correction = p.adjust(p_value, method = correction)) %>%
  arrange(correction) %>%
  rename(p_adjusted = `correction`)
print(dmin, n=10)

# Filter by significant excess allele sharing
dmin_sign <- dmin %>%
  filter(p_adjusted < 0.05)
print(dmin_sign, n=nrow(dmin_sign))

write_delim(dmin_sign, paste("Dmin_sign_", dataset, ".csv", sep = ""), delim = "\t")

## D-stats for BBAA topology
bbaa <- read_tsv(paste("hyp_", dataset, "_dtrios_BBAA.txt", sep ="")) %>%
  rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
  mutate(correction = p.adjust(p_value, method = correction)) %>%
  arrange(correction) %>%
  rename(p_adjusted = `correction`)
print(bbaa, n=20)

write_delim(bbaa, paste("BBAA_", dataset, ".csv", sep = ""), delim = "\t")

# Filter by significant excess allele sharing
bbaa_sign <- bbaa %>%
  filter(p_adjusted < 0.05)
print(bbaa_sign, n=nrow(bbaa_sign))

write_delim(bbaa_sign, paste("BBAA_sign_", dataset, ".csv", sep = ""), delim = "\t")

## Overlap Dmin, BBAA
semi_join(dmin_sign[1:3], bbaa_sign[1:3])


## Pairs with evidence for gene flow
bbaa_pairs <- (unique(bbaa_sign[2:3]))
print(bbaa_pairs, n=nrow(bbaa_pairs))

dmin_pairs <- (unique(dmin_sign[2:3]))
print(dmin_pairs, n=nrow(dmin_pairs))


## Proportions
species <- readLines("species_order_alpha.txt")[2:15]
trios <- combinations(n = length(species), r = 3, v = species)
duos <- combinations(n = length(species), r = 2, v = species)