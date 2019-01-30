#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla sample_assignment_msmc.R distribute_samples_msmc_and_cc.R cross_cc.R ~/Desktop/chapter2/metadata/sample_info.txt ./msmc
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c('msmc/distribute_samples_msmc_and_cc.R', 'msmc/cross_cc.R','~/Desktop/chapter2/metadata/sample_info.txt','msmc/msmc')
args = commandArgs(trailingOnly=FALSE)
args = args[7:10]
print(args)
# setup -----------------------

library(tidyverse)

distr_script <- as.character(args[1])
cross_script <- as.character(args[2])
sample_info <- as.character(args[3])
prefix <- as.character(args[4])
source(distr_script)
source(cross_script)
# -----------------------------------
data <- read_delim(sample_info, delim = '\t') %>% 
  select(label:geo) %>% 
  filter( !(spec %in% c('flo','tab','tor')))
  
group_layout <- tibble(n = 10:13,
                       grps_msmc_n = ceiling(n/4),
                       n_4 = c(1:3,1),
                       n_3 = c(2,1,0,3)#, check = 4*four+3*three, check_n = four+three
  )
  
group_sizes <- data  %>%
  group_by(spec,geo) %>%
  summarise( n = length(geo),
             inds = list(label)) %>%
  left_join(.,group_layout) %>% 
  ungroup()

set.seed(27678)
msmc_grouping <- group_sizes %>%
  select(inds,n_4,n_3,spec,geo) %>%
  pmap(collapse_samples_msmcs) %>% 
  bind_rows() %>%
  mutate(msmc_run = row_number()) %>%
  select(msmc_run,spec:group_nr,group_size,samples)

cc_grouping <- group_sizes %>%
  select(inds,spec,geo) %>%
  pmap(collapse_samples_cc) %>% 
  bind_rows() %>%
  select(spec:group_nr,samples)

cc_samples <- cc_grouping %>%
  mutate(id = str_c(spec,'_',geo,'_',group_nr)) %>% 
  select(id,samples)

cc_tibbles <- cc_grouping %>%
  group_by(spec,geo) %>%
  count() %>%
  ungroup() %>%
  group_by(geo) %>%
  summarise(content = list(tibble(spec=spec,n=n))) %>%
  bind_cols(., .$content %>% map(paste_groups_cc) %>% tibble) %>%
  setNames(.,nm = c('geo','content', 'contrasts'))
 
cc_output <- cc_tibbles %>% 
  select(-content) %>% 
  pmap(spread_tibble_cc) %>% 
  bind_rows() %>%
  mutate( id1 = str_c(spec_1,'_',geo,'_',group_1),
          id2 = str_c(spec_2,'_',geo,'_',group_2)) %>%
  left_join(.,cc_samples %>% set_names(.,nm = c('id1','samples_1')))%>%
  left_join(.,cc_samples %>% set_names(.,nm = c('id2','samples_2'))) %>% 
  mutate(run_nr = row_number()) %>%
  select(run_nr,geo,spec_1,spec_2,contrast_nr,samples_1,samples_2)

write_delim(x = msmc_grouping,path = str_c(prefix,'_grouping.txt'),delim = '\t')
write_delim(x = cc_output,path = str_c(prefix,'_cc_grouping.txt'),delim = '\t')
