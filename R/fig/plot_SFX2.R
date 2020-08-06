#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFX2.R 2_analysis/ccf/
# ===============================================================
# This script produces Suppl. Figure X of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/ccf/" )
# script_name <- "R/fig/plot_SFX.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name)))
args = args[7:length(args)]
cat(' ')
cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args),': ',crayon::green(args),'\n'))
cli::rule(right = getwd())

# config -----------------------
ccf_dir <- as.character(args[1])

ccf_files <- dir(ccf_dir, pattern = "ccf.LG20.18180puebel")
  
read_ccf <- function(path = ccf_dir, file = ccf_files[[1]]){
  ind1 <- file %>% str_remove("^ccf.LG[0-9]{2}.") %>% str_split(pattern = "\\.") %>% .[[1]] %>% .[[1]]
  ind2 <- file %>% str_remove("^ccf.LG[0-9]{2}.") %>% str_split(pattern = "\\.") %>% .[[1]] %>% .[[2]]
  vroom::vroom(str_c(path,file)) %>%
    set_names(nm = c("CHROM", "POS", "PostMedian",
                     str_c(ind1, "_1-",ind1, "_2"),
                     str_c(ind1, "_1-",ind2, "_1"),
                     str_c(ind1, "_1-",ind2, "_2"),
                     str_c(ind1, "_2-",ind2, "_1"),
                     str_c(ind1, "_2-",ind2, "_2"),
                     str_c(ind2, "_2-",ind2, "_2"))) %>%
    pivot_longer(names_to = "pair", 
                 cols = -(CHROM:PostMedian)) %>%
    mutate(prep = pair) %>%
    separate(prep, into = c("ind1","ind2"), sep = "-") %>%
    separate(ind1, into = c("ind1","hap1"), sep = "_",convert = TRUE)%>%
    separate(ind2, into = c("ind2","hap2"), sep = "_",convert = TRUE) %>%
    mutate(spec1 = ind1 %>% str_sub(-6,-4),
            spec2 = ind2 %>% str_sub(-6,-4),
           value = if_else(hap1 == 1, value, -value)) %>%
    arrange(pair, PostMedian) %>%
    group_by(pair) %>%
    mutate(check = lag(value, default = 0) == value & lead(value, default = 0) == value) %>%
    ungroup() %>%
    filter(!check)
}

ccf_data <- cff_files %>% 
  purrr::map(read_ccf, path = ccf_dir) 

ccf_data  %>%
  ggplot(aes(x = PostMedian, y = value, group = pair, color = spec2)) +
  geom_line()+
  facet_grid(hap1 ~ .)+
  scale_x_log10()+
  scale_color_manual(values = clr)
