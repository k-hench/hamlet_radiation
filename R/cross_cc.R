sample_groups_cc <- function(size,max_n){sample(x = 1:size,size = max_n,replace = FALSE)}
sample_groups_cc_vector <- function(size,max_n){
  purrr::pmap(tibble(size=size,max_n = max_n),sample_groups_cc) 
}

cross_groups_cc <- function(input_tibble){
  cross_df(.l = list(pop1 = input_tibble$spec, pop2 = input_tibble$spec),.filter = `>=`) %>%
    left_join(input_tibble %>% setNames(.,nm = c('pop1','size1'))) %>%
    left_join(input_tibble %>% setNames(.,nm = c('pop2','size2'))) %>%
    mutate(max_n = ifelse(size1 > size2, size2, size1),
           groups_1 = sample_groups_cc_vector(size = size1,max_n = max_n),
           groups_2 = sample_groups_cc_vector(size = size2,max_n = max_n))  %>% 
    select(-size1,-size2)
}

extract_groups_cc <- function(pop1,pop2,groups_1,groups_2){
  str_c(str_c(pop1,groups_1,sep = '.'), 
        str_c(pop2,groups_2,sep = '.'),
        sep = '-')#,collapse = ', ')
}

paste_groups_cc <- function(input_tibble){
  input_tibble %>%
    cross_groups_cc %>% 
    select(-max_n) %>% 
    pmap(.,extract_groups_cc)  %>% 
    unlist() 
}

spread_tibble_cc <- function(geo,contrasts){
  tibble(contrast = contrasts) %>%
    mutate(geo = geo,
           contrast_nr = row_number(),
           pre = contrast) %>%
    separate(pre,into = c('group_1', 'group_2'),sep = '-') %>%
    separate(group_1,into = c('spec_1', 'group_1'),sep = '\\.') %>%
    separate(group_2,into = c('spec_2', 'group_2'),sep = '\\.')
}
