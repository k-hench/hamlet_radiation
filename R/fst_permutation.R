library(tidyverse)
library(vroom)
library(ggtext)

rand_path <- "2_analysis/fst_signif/random/"
rand_files <- dir(path = rand_path)

get_random_fst <- function(file){
  nm <- file %>% str_remove(pattern = ".*\\/") %>% str_remove("_random_fst.tsv.gz")
  vroom::vroom(file = file,
               delim = "\t",
               col_names = c("idx", "type", "mean_fst", "weighted_fst"),
               col_types = "dcdd") %>%
  mutate(run = nm)
}

data <- str_c(rand_path, rand_files) %>%
  map_dfr(.f = get_random_fst)

get_percentile <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  
  sprintf("%.3f", sum(ran < real_fst) / length(ran))
}

get_n_above<- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  
  sum(ran > real_fst)
}

get_n_total <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
 length(ran)
}

data_grouped <- data %>% 
  group_by(run) %>%
  nest() %>%
  ungroup() %>%
  mutate(percentile = map_chr(data, get_percentile),
         above = map_dbl(data, get_n_above),
         total = map_dbl(data, get_n_total),
         label = str_c(percentile, "<br>(",above, "/", total, ")"))

data %>%
  filter(type == "random") %>%
  ggplot() +
  geom_vline(data = data %>%
               filter(type == "real_pop"),
             aes(xintercept = weighted_fst),
             color = "red") +
  geom_density(aes( x = weighted_fst ),
               fill = rgb(0,0,0,.3)) +
  geom_richtext(data = data_grouped,
                aes(x = .08, y = 400, label = label), 
                hjust = .5, vjust = 1, size = 3,
                label.padding = unit(3,"pt"),
                label.size = 0,
                label.color = "transparent",
                fill = "transparent") +
  facet_wrap(run ~ .) +
  theme_minimal()

scl <- .75
ggsave("~/Desktop/fst_permutation.pdf",
       width = 16 * scl,
       height = 9 * scl,
       device = cairo_pdf)
