take_inds <- function(inds,n){
  new_sample <- sample(x = inds$remain, size = n, replace = FALSE)
  inds_remain <- dplyr::setdiff(inds$remain, new_sample)
  list(remain = inds_remain, 
       inds_sample = rlist::list.append(inds$inds_sample,new_sample))
}

sample_size_class <- function(input_inds, n_samples, size_samples){
  while (length(input_inds$remain) >= size_samples & n_samples > 0) {
    n_samples <- n_samples -1
    input_inds <- take_inds(input_inds,size_samples)
  }
  input_inds
}

compose_samples_msmcs <- function(inds,n_4,n_3){
  ind_list <- list(remain = inds,inds_sample=list())
  sample_size_class(ind_list,n_samples = n_4,size_samples = 4) %>%
    sample_size_class(.,n_samples = n_3,size_samples = 3)
}

compose_samples_cc <- function(inds){
  n_samples <- length(inds)%/%2
  ind_list <- list(remain = inds,inds_sample=list())
  sample_size_class(ind_list,n_samples = n_samples,size_samples = 2)
}

collapse_samples_msmcs <- function(inds,n_4,n_3,spec,geo){
  compose_samples_msmcs(inds = inds, n_4 = n_4, n_3 = n_3)$inds_sample %>% 
    map(.,str_c,collapse=', ') %>% 
    unlist() %>% 
    tibble(samples = .) %>%
    mutate(group_size = str_count(samples,',')+1,
           spec = spec,
           geo = geo,
           group_nr = row_number())
}

collapse_samples_cc <- function(inds,spec,geo){
  compose_samples_cc(inds = inds)$inds_sample %>% 
    map(.,str_c,collapse=', ') %>% 
    unlist() %>% 
    tibble(samples = .) %>%
    mutate(spec = spec,
           geo = geo,
           group_nr = row_number())
}