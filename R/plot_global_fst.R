#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla plot_global_fst.R ${globals.file.txt}
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c('hon.logs.txt','fst_functions.R','~/Desktop/chapter2/R/project_config.R')
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]
print(args)
library(tidyverse)
library(hypoimg)
globals_file <- as.character(args[1])
functions_script <- as.character(args[2])
proj_config <- as.character(args[3])
source(functions_script)
source(proj_config)

globals <- read_delim(globals_file, delim = '\t',
                      col_names = c('loc','run','mean','weighted')) %>%
  mutate(loc_run = str_c(loc,'_',run),
         run = fct_reorder(run,weighted),
         loc_run = fct_reorder(loc_run,weighted)) 

plot_hypo <- function(loc_run,left,right,circle_fill,...){ 
  tibble(loc_run = loc_run,
         grob = list( hypo_anno_pair(left = left,
                                     right = right, 
                                     circle_fill = circle_fill,
                                     ...) %>% ggplot2::ggplotGrob()
                      )
         )
  }

hypo_tab <- globals %>% 
  mutate(pre = run, circle_fill = clr_loc[loc]) %>% 
  separate(pre, into=c('left','right')) %>%
  mutate(loc_run = as.character(loc_run),
         left = sp_names[left], 
         right = sp_names[right]) %>%
  select(loc_run,left:circle_fill) %>%
  purrr::pmap(plot_hypo,circle_color = 'black',circle_lwd = .2) %>%
  bind_rows()

hypo_tab$loc_run <- refactor(hypo_tab,globals)

locs <- tibble(geo = globals$loc %>% as.factor() %>% levels(),
               color_map = clr_loc[geo])
               
legend_flag <- hypo_legend_flag_single(geo = loc_names[locs$geo] %>% tolower(),
                      color_map = locs$color_map,
                      plot_names = TRUE,
                      circle_color = 'black')
hypo_scale <- 2
p1 <- globals %>% 
  left_join(hypo_tab) %>%
  mutate(x = as.numeric(loc_run)) %>%
  ggplot(aes(x = .5, y = weighted)) +
  facet_grid(.~loc_run)+
  geom_bar(stat = 'identity', aes(fill = loc)) +
  geom_hypo_grob( aes(x = .5,grob = grob),y = .05,angle=-90,height = hypo_scale,width = hypo_scale)+
  scale_y_continuous(name = expression(Genome~wide~weighted~italic('F'[ST])),
                     expand = c(0,0),
                     limits = c(-max(globals$weighted)/8,max(globals$weighted)))+
  scale_x_discrete(labels = globals$run,expand=c(.05,.05))+
  scale_fill_manual(values = clr_loc,guide=FALSE)+
  theme_bw(base_size = 10, base_family = "Helvetica") %+replace% 
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),panel.spacing = unit(2,'pt'), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90),
        #axis.line = element_line(),
        #axis.line.y = element_line(color = hypogen::hypo_clr_lg),
        axis.line.y.right = element_line(color = hypogen::hypo_clr_lg),
        strip.background = element_rect(fill = NA,color = hypogen::hypo_clr_lg),
        legend.background = element_rect(fill = "transparent", 
                                         color = NA),
        legend.key = element_rect(fill = "transparent",color = NA),
        strip.text = element_blank(),
        legend.position = c(.1,.9))

plot_complete <- ggdraw(p1)+
  draw_plot(legend_flag,.0125,.7,.3,.3)

ggsave(plot = plot_complete,filename = 'global_fst.pdf',width = 8,height = 9,device = cairo_pdf)
