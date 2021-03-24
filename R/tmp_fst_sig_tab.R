read_tsv("~/Downloads/fst_fdr_by.tsv") %>% 
  select(run, p_perm, loc_n, fdr_correction_factor, fdr_alpha ) %>% 
  mutate(`050` = .05/ fdr_correction_factor,
         `010` = .01/ fdr_correction_factor,
         `001` = .001/ fdr_correction_factor,
         sig = ifelse(p_perm < `001`, "***", ifelse(p_perm < `010`, "**", ifelse(p_perm < `050`, "*", "-")))) %>% 
  filter(grepl( "hon", run)) 
