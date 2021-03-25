#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFY.R 2_analysis/fst_signif/random/
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/fst_signif/random/")
# script_name <- "R/fig/plot_F1.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(vroom)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
rand_path <- as.character(args[1])
rand_files <- dir(path = rand_path)

get_random_fst <- function(file){
  nm <- file %>% str_remove(pattern = ".*\\/") %>%
    str_remove("_random_fst.tsv.gz")
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
  
  sprintf("%.7f", sum(ran < real_fst) / length(ran))
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
  mutate(real_pop = map_dbl(data, function(data){data$weighted_fst[data$type == "real_pop"]}),
         percentile = map_chr(data, get_percentile),
         above = map_dbl(data, get_n_above),
         total = map_dbl(data, get_n_total),
         label = str_c(sprintf("%.2f", as.numeric(percentile)), "<br>(",above, "/10<sup>",log10(total),"</sup>)")) %>%
  arrange(-as.numeric(percentile), -real_pop) %>%
  mutate(rank = row_number(),
         run = fct_reorder(run, rank)) 

p_done <- data %>%
  mutate(run = factor(run, levels = levels(data_grouped$run))) %>%
  filter(type == "random") %>%
  ggplot() +
  geom_vline(data = data %>%
               mutate(run = factor(run, levels = levels(data_grouped$run))) %>%
               filter(type == "real_pop"),
             aes(xintercept = weighted_fst),
             color = "red") +
  geom_density(aes( x = weighted_fst ),
               fill = rgb(0,0,0,.3)) +
  geom_richtext(data = data_grouped,
                aes(x = .08, y = 470, label = label), 
                hjust = .5, vjust = 1, size = 3,
                label.padding = unit(3,"pt"),
                label.size = 0,
                label.color = "transparent",
                fill = "transparent") +
  scale_x_continuous(breaks = c(0, .05, .1)) +
  facet_wrap(run ~ ., ncol = 4,dir = "v") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

scl <- 1.4
hypo_save(p_done, filename = 'figures/SFY.pdf',
          width = f_width * scl,
          height = f_width * .6  * scl,
          device = cairo_pdf,
          comment = plot_comment,
          bg = "transparent")

  data_export <- data_grouped %>%
  select(run, real_pop, percentile) %>%
  mutate(pre = run,
         p_perm = 1 - as.numeric(percentile)) %>%
  separate(pre, into = c("p2", "p1"))

write_tsv(x = data_export, file = "2_analysis/summaries/fst_permutation_summary.tsv")