#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla plot_fst.R ${location} ${globals.file.txt} fst_functions.R project_config.R
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c('hon','hon.logs.txt','fst_functions.R','~/Desktop/chapter2/R/project_config.R')
args = commandArgs(trailingOnly=FALSE)
args = args[7:10]
print(args)
# setup -----------------------
library(tidyverse)
library(hypogen)
library(hypoimg)

# config -----------------------
loc <- as.character(args[1])
globals_file <- as.character(args[2])
functions_script <- as.character(args[3])
proj_config <- as.character(args[4])
source(functions_script)
source(proj_config)
# load data -------------------
files <- dir(pattern = '.50k.windowed.weir.fst.gz')
run_files <- files %>% str_remove(str_c(loc,'-')) %>% str_remove(str_c('.50.*'))

globals <- read_delim(globals_file, delim = '\t',
                      col_names = c('loc','run','mean','weighted')) %>%
  mutate(loc_run = str_c(loc,'_',run),
         run = fct_reorder(run,weighted),
         loc_run = fct_reorder(loc_run,weighted))

data <- purrr::pmap(tibble(file = files,run = run_files),
                    hypo_import_windows) %>%
  bind_rows() %>%
  mutate(loc_run = str_c(loc,'_',RUN))

# global bars ------------------------
global_bar <- globals %>%
  select(weighted,loc_run) %>%
  mutate(loc_run = as.character(loc_run)) %>%
  setNames(.,nm = c('fst','loc_run')) %>%
  pmap(.,fst_bar_row) %>%
  bind_rows()

# annotaton ---------------------------
runs <- globals %>%
  group_by(loc_run) %>%
  count() %>% ungroup() %>%
  select(-n) %>%
  mutate(pre = loc_run) %>%
#  mutate(loc_run = str_remove(as.character(loc_run),'^..._')) %>%
  separate(pre,into = c('loc','right_short','left_short')) %>%
  mutate(loc_run = as.character(loc_run),
    left = sp_names[left_short],
         right = sp_names[right_short],
         circle_fill_left = clr[left_short],
         circle_fill_right = clr[right_short])

grob_tibble <- runs %>% select(1:2,5:8) %>% pmap(.,plot_pair) %>% bind_rows()

# arrange factor order ----------------
data$loc_run <- refactor(data, globals)
grob_tibble$loc_run <- refactor(grob_tibble, globals)
global_bar$loc_run <- refactor(global_bar,globals)

# plotting ---------------
p1_1 <- ggplot()+
  facet_grid( loc_run~., as.table = TRUE)+
  geom_hypo_LG()+
  geom_point(data = data, aes(x = GPOS, y = WEIGHTED_FST),
	           size = plot_size, color = plot_clr) +
  scale_fill_hypo_LG_bg() +
  scale_x_hypo_LG(name = loc_names[loc])+
  scale_y_continuous(name = expression(italic('F'[ST])),limits = c(-.1,1))+
 # scale_color_manual(values = clr)+
  theme_bw(base_size = 10) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        strip.background = element_rect(fill = NA,color = hypo_clr_lg),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent",color = NA),
        strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text())

x_fst <- seq(0,ceiling2d(max(globals$weighted)),length.out = 5)[1:4]
p1_2 <- ggplot()+
  facet_grid( loc_run~., as.table = TRUE)+
  geom_hypo_grob(data = grob_tibble,
                 aes(grob = grob, x = .5,y = .6),
                 angle = 0, height = .5)+
  geom_rect(data = global_bar,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),col=NA , fill = plot_clr)+
  scale_x_continuous(name = expression(global~italic('F'[ST])),
                     limits = c(0,1),expand = c(0,0),
                     breaks = rescale_fst(x_fst),labels = x_fst,position = "top")+
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "y2"),limits = c(-.1,1))+
 # scale_color_manual(values = clr)+
  theme_bw(base_size = 10) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.y = element_line(color = hypogen::hypo_clr_lg),
        axis.line.y.right = element_line(color = hypogen::hypo_clr_lg),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA,color = hypo_clr_lg),
        legend.background = element_rect(fill = "transparent",color = NA),
        legend.key = element_rect(fill = "transparent",color = NA),
        strip.text = element_blank(),
        legend.position = 'none')

p1 <- cowplot::plot_grid(p1_1,p1_2,ncol = 2,align = 'h',rel_widths = c(1,.15))
n_runs <- runs$loc_run %>% length()
ggsave(p1, filename = str_c('fst_',loc,'.pdf'),width = 297*.95,height = (15+20*n_runs)*.95,units = 'mm')
