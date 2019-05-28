library(tidyverse)

read_tsv('pop_prep.tsv', col_names = c('ID','POP')) %>%
  mutate(nr = POP %>% factor() %>% as.numeric(),
    group = sample(c('A','B')[nr], replace = FALSE)) %>%
  select(ID,group) %>%
  write_tsv('random_pop.txt',col_names = FALSE)