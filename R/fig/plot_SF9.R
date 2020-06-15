#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig_rev/plot_SF9.R 2_analysis/GxP/50000/
# ===============================================================
# This script produces Suppl. Figure 9 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/GxP/50000/')
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
gxp_path <- as.character(args[1])
# load data -------------------
trait_tib  <- tibble(file = dir(gxp_path) %>% .[str_detect(.,"Bars|Peduncle|Snout")]) %>%
  mutate(prep = file) %>%
  separate(prep , into = c("trait", "model_type", "win", "step", "filetype", "zip"),
           sep = "\\.") %>%
  select(file, trait, model_type) %>%
  mutate(path = gxp_path)

data <- pmap_dfr(trait_tib,get_gxp_both_models)

p <- data %>%
  ggplot(aes(x = gpos, y = AVG_p_wald))+
  geom_hypo_LG()+
  geom_point(color = plot_clr, size = .3)+
  scale_x_hypo_LG()+
  scale_fill_hypo_LG_bg()+
  labs(y = expression(G~x~P~(average~italic(p)[wald])))+
  facet_grid(trait+model_type ~ ., scales = "free_y")+
  theme_hypo()

#scl <- 1
ggsave(filename = "figures/SF9.png",
       plot = p,
       width = 11,
       height = 7,
       dpi = 600
       )
