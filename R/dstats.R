#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla dstats.R

### Extract significant Dmin, BBAA trios
### -------------------------------------

### Jul 12, 2021, Martin Helmkampf

library(tidyverse)

dataset <- "ld05"
# Alternatives: 5kb, ld06, ld05

correction <- "holm"
# Alternatives: bonferroni, fdr


## D-stats for topologies minimizing D (Dmin)
dmin <- read_tsv(paste("hyp_", dataset, "_dtrios_Dmin.txt", sep ="")) %>%
  rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
  mutate(correction = p.adjust(p_value, method = correction)) %>%
  arrange(correction) %>%
  rename(!!correction := correction)
print(dmin, n=10)

# Filter by significant excess allele sharing
dmin_sign <- dmin %>%
  filter(holm < 0.05)
print(dmin_sign, n=nrow(dmin_sign))

write_delim(dmin_sign, paste("Dmin_sign_", dataset, ".csv", sep = ""), delim = "\t")

## D-stats for BBAA topology
bbaa <- read_tsv(paste("hyp_", dataset, "_dtrios_BBAA.txt", sep ="")) %>%
  rename(Z_score = `Z-score`, p_value = `p-value`, f4_ratio = `f4-ratio`) %>%
  mutate(correction = p.adjust(p_value, method = correction)) %>%
  arrange(correction) %>%
  rename(!!correction := correction)
print(bbaa, n=20)

# Filter by significant excess allele sharing
bbaa_sign <- bbaa %>%
  filter(holm < 0.05)
print(bbaa_sign, n=nrow(bbaa_sign))

write_delim(bbaa_sign, paste("BBAA_sign_", dataset, ".csv", sep = ""), delim = "\t")


## Overlap Dmin, BBAA
semi_join(dmin_sign[1:3], bbaa_sign[1:3])


## Pairs with evidence for gene flow (BBAA)
pairs <- (unique(bbaa_sign[2:3]))
print(pairs, n=nrow(pairs))
