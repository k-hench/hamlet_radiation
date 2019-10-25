library(tidyverse)
custom_pkg <- c("GenomicOriginsScripts", "hypoimg", "hypogen")
all_pkg <- system("grep -r ^library R/* | sed 's/.*library(//g; s/)//' | sort | uniq", intern = TRUE) 


pkg_table <- function(pkg){
  tibble(pkg = pkg,
         local_version = packageVersion(pkg)%>% str_c())
}


pkg_citation <- function(pkg, width = 50){
  citation(package = pkg) %>% 
             toBibtex() %>%
             str_replace(., pattern = "(@.*)\\{,", replacement = str_c('\\1\\{',pkg,',')) %>%
             str_wrap(., width, 2, 4) %>%
             str_c(.,collapse = '\n')
}

cluster_pkg <- tibble(pkg = c('FastEPRR', 'logisticPCA'),
                      cluster_version = c('1.0', '0.2.9000'))

all_pkg[! (all_pkg %in% custom_pkg)] %>%
  c(custom_pkg,.) %>%
  map_dfr(pkg_table) %>%
  arrange(pkg) %>%
  left_join(cluster_pkg) %>%
  mutate(cite_version = if_else(is.na(cluster_version), 
                           local_version, 
                           cluster_version),
#         output = str_c('\\textsw{', pkg, '} \\citef[', pkg, ']{} (', cite_version, ')')
         output = str_c('\\textsw{', pkg, '} (', cite_version, ')')
         ) %>%
  select(output) %>% 
  unname() %>% 
  unlist() %>% 
  str_c(collapse = ', ') %>% 
  cat()
  

all_pkg[! (all_pkg %in% custom_pkg)] %>%
  c(custom_pkg,.) %>%
  map(pkg_citation) %>%
  str_c(collapse = '\n\n') %>% 
  cat()

