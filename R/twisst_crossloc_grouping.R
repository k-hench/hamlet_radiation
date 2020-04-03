#!/usr/bin/env Rscript
# Rscript --vanilla twisst_crossloc_grouping.R
library(tidyverse)

geos <- c("bel", "hon", "pan")
specs <- c("nig", "pue", "uni")

filter <- function(x, y) x >= y

geo_cross <- cross_df(list(geo1 = geos, geo2 = geos),.filter = filter) %>%
  unite("geo_cross", geo1,geo2)

spec_cross <- cross_df(list(spec1 = specs, spec2 = specs),.filter = filter) %>%
  unite("spec_cross", spec1,spec2)

full_cross <- cross_df(list(geo_cross = geo_cross$geo_cross, spec_cross = spec_cross$spec_cross)) %>%
  separate(geo_cross, into = c("geo1", "geo2"))%>%
  separate(spec_cross, into = c("spec1", "spec2")) %>%
  mutate(runnr = row_number(),
         pop1 = str_c(spec1, geo1),
         pop2 = str_c(spec2, geo1),
         pop3 = str_c(spec1, geo2),
         pop4 = str_c(spec2, geo2))

full_cross %>%
  select(runnr,
         starts_with("pop")) %>%
  write_tsv(path = "twisst_runs_config.tsv")