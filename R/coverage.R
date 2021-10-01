library(tidyverse)
library(GenomicOriginsScripts)
library(prismatic)

cov_path <- "1_genotyping/0_dedup_bams/cov/"
files <- dir(cov_path)

import <- function(file, cov_path =  "1_genotyping/0_dedup_bams/cov/"){
  file_path <- str_c(cov_path, file) 

  vroom::vroom(file_path) %>% 
    mutate(file = file)
}

data <- files %>%
  map_dfr(import) %>%
  mutate(depth_w = covbases * meandepth)

data %>% 
  group_by(file) %>% 
  summarise(depth = sum(depth_w) / sum(covbases)) %>% 
  mutate( spec = file %>% str_remove(pattern = ".both") %>%
            str_remove(pattern = ".cov") %>% 
            str_remove(pattern = ".dedup") %>% 
            str_remove(pattern = "\\.[0-9]") %>% 
            str_sub(start = -6, end = -4),
          loc = file %>% str_remove(pattern = ".both") %>%
            str_remove(pattern = ".cov") %>% 
            str_remove(pattern = ".dedup") %>% 
            str_remove(pattern = "\\.[0-9]") %>%  str_sub(start = -3, end = -1),
          grp = str_c(spec, loc),
          fll = forcats::fct_reorder(file %>%
                                       str_remove(pattern = ".cov") %>% 
                                       str_remove(pattern = ".dedup") %>% 
                                       str_remove(pattern = "\\.[0-9]") %>% str_remove("[a-z]{6}"), -depth)) %>% 
  dplyr::select(-fll) %>% 
  arrange(grp, -depth) %>% write_tsv("~/Desktop/coverage.tsv")
  ggplot(aes(fll, y = depth, fill = spec)) +
  facet_wrap(grp ~ ., scales = "free") +
  geom_bar(stat = "identity", aes(color = after_scale(clr_darken(fill)))) +
  # coord_cartesian(ylim = c(95, 101)) +
  scale_fill_manual(values = clr) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1))

ggsave("~/Desktop/coverage.pdf", width = 9, height = 10, device = cairo_pdf)
