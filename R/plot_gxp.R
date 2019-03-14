#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla plot_gxp.R $BASE_DIR/R/project_config.R
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c('~/Desktop/chapter2/R/project_config.R')
args = commandArgs(trailingOnly=FALSE)
args = args[7]
print(args)
# setup -----------------------
library(tidyverse)
library(hypogen)
library(hypoimg)

# config -----------------------
proj_config <- as.character(args[1])
source(proj_config)
# load data -------------------
files <- dir(pattern = '.lm.50k.5k.txt.gz')
traits <- files %>%  str_remove(str_c('.lm.*'))

data <- purrr::pmap(tibble(file = files,run = traits),
                    hypo_import_windows) %>% 
  bind_rows() %>%
  group_by(RUN) %>%
  mutate(ord = max(AVG_p_wald)) %>%
  ungroup() %>%
  mutate(RUN = fct_reorder(RUN,ord)) %>% 
  select(-ord)


fish <- tibble(spec = sp_names[traits[!traits %in%  c(hypo_trait_img$trait,'PC_d1')]],
       RUN = traits[!traits %in% c(hypo_trait_img$trait,'PC_d1')]) %>%
  left_join(.,hypo_img) %>% mutate(grob = l)%>% 
  select(RUN,grob)

img_grobs <- hypoimg::hypo_trait_img %>% 
  mutate(trait = ifelse(trait=='PC_1d','PC_d1',trait),
         RUN = trait,
         grob = grob_circle) %>% 
  select(RUN,grob) %>%
  bind_rows(., fish) %>% 
  mutate(RUN = factor(RUN,levels = levels(data$RUN)))

p1 <- ggplot()+
  geom_hypo_LG()+
  geom_point(data = data, aes(x = GPOS, y = AVG_p_wald),
             size = plot_size, color = plot_clr) +
  geom_hypo_grob(data = img_grobs,
                 aes(grob = grob, x = .965, y=.5, 
                     angle = 0, height = .85))+
  scale_fill_hypo_LG_bg() +
  scale_x_hypo_LG(limits = c(0,5.95e+08))+
  facet_grid( RUN~., as.table = TRUE, scales = 'free_y')+
  scale_y_continuous(name = expression(-log[10](italic(p)~value)))+
  theme_hypo()+
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(plot = p1, filename = 'gxp.pdf',width = 297*.95,height = (10+20*length(traits))*.95,units = 'mm')


