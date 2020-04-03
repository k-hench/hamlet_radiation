#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF1.R 2_analysis/summaries/fst_globals.txt \
# ressources/other_studies/stankowski_etal_2019.tsv \
# ressources/other_studies/han_etal_2017.tsv
# ===============================================================
# This script produces Figure 1 of the study "The genomic onset of a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/summaries/fst_globals.txt',
#           "ressources/other_studies/stankowski_etal_2019.tsv",
#           "ressources/other_studies/han_etal_2017.tsv")
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(tidyverse)
library(ggstance)
library(GenomicOriginsScripts)
library(scales)
library(hypoimg)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
globals_file <- as.character(args[1])
monkeyflower_file <- as.character(args[2])
finches_file <- as.character(args[3])
# load data -------------------
fst_range <- c(87.690, 208.643)

heliconius <- tibble(run = c("cyd-pac", "cyd-mel", "pac-mel"),
       fst_y = c(96.620, 121.756, 131.725)) %>%
  mutate(mean_fst = rescale(fst_y, to = c(0,.8), from = fst_range),
         study = "*Heliconius*<sup>1</sup>") 

sunflowers <- tibble(run = c("ann-pet", "ann-deb", "deb-arg", "pet-arg"),
                     mean_fst = c(.3, .35, .51, .48),
                     study = "Sunflowers<sup>2</sup>")

hamlets <- read_tsv(globals_file,
                    col_names = c('loc','run','mean_fst','weighted_fst'))%>%
  mutate(study = "Hamlets")

monkeyflowers <- read_tsv(monkeyflower_file,
                          col_names = c('run','mean_fst')) %>%
  mutate(study = "Monkeyflowers")
  
finches <- read_tsv(finches_file) %>%
  select(`Species pairs`, `Genome-wide mean F ST`) %>%
  set_names(nm = c("run", "mean_fst")) %>%
  mutate(study = "Finches")

flycatchers <- tibble(run = c("close", "far","alb-hyp", "alb-spe", "alb-sem", "hyp-spe", "hyp-sem", "spe-sem"),
                      mean_fst = c(0.023, 0.041, 0.274, 0.201, 0.326, 0.322, 0.398, 0.394),
                      study = "Flycatchers")

data <- hamlets %>% select(run, mean_fst, study) %>%
  bind_rows(monkeyflowers %>% select(run, mean_fst, study)) %>%
  bind_rows(finches) %>%
  bind_rows(heliconius) %>%
  bind_rows(sunflowers) %>%
  bind_rows(flycatchers) %>%
  mutate(study = fct_reorder(study,.x = mean_fst, .fun = median,.desc = TRUE))

data_sumary <- data %>% 
  group_by(study) %>% 
  summarise(n = length(mean_fst),
            lower = quantile(mean_fst, .25),
            upper = quantile(mean_fst, .75),
            min = min(mean_fst),
            max = max(mean_fst))

bin_wdh <- .35

img_hel <- hypo_read_svg("ressources/img/heliconius.c.svg")
img_sun <- hypo_read_svg("ressources/img/sunflower.c.svg")
img_fin <- hypo_read_svg("ressources/img/finch.c.svg")
img_fly <- hypo_read_svg("ressources/img/flycatcher.c.svg")
img_mon <- hypo_read_svg("ressources/img/monkeyflower.c.svg")

img_hyp <- (hypoimg::hypo_outline %>%
              ggplot(aes(x,y))+
              geom_polygon(fill = "black")+
              coord_equal()+
              theme_void() ) %>% ggplotGrob()
p <- data %>%
  ggplot(aes(y = study))+
  annotation_custom(grob = img_hyp, xmax = 1.35, xmin = 1, ymin = 5.55, ymax = 6.45)+
  annotation_custom(grob = img_hel, xmax = 1.35, xmin = 1, ymin = 4.5, ymax = 5.5)+
  annotation_custom(grob = img_fin, xmax = 1.35, xmin = 1, ymin = 3.55, ymax = 4.45)+
  annotation_custom(grob = img_fly, xmax = 1.35, xmin = 1, ymin = 2.55, ymax = 3.45)+
  annotation_custom(grob = img_mon, xmax = 1.35, xmin = 1, ymin = 1.5, ymax = 2.5)+
  annotation_custom(grob = img_sun, xmax = 1.35, xmin = 1, ymin = .5, ymax = 1.5)+
  geom_boxploth(aes(x = mean_fst, y = study),
                width = bin_wdh, fill = clr_below) +
  geom_text(data = data_sumary, 
            aes(y = study, x = .825, label = str_c("italic(n):~",n)),
            parse = TRUE, hjust = 0)+
  scale_x_continuous(name = "Genome-wide differentiation (mean *F<sub>ST</sub>*)",
                     limits = c(0,1.3))+
  theme_minimal()+
  labs(caption = "1: median of 5kb windows on autosomes\n2: based on transcriptomes")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = clr_below, size =.2),
        plot.caption = element_text(color = rgb(.4, .4, .4)))

hypo_save(filename = 'figures/SF1.pdf',
          plot = p,
          width = 6,
          height = 3,
          device = cairo_pdf,
          comment = plot_comment)