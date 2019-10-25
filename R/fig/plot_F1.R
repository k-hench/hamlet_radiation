#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F1.R 2_analysis/dxy/50k/ 2_analysis/fst/50k/
# ===============================================================
# This script produces Figure 1 of the study "The genomic origins of a marine radiation"
# by Hench, McMillan an Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/dxy/50k/', '2_analysis/fst/50k/')
args <- commandArgs(trailingOnly=FALSE)
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
dxy_dir <- as.character(args[1])
fst_dir <- as.character(args[2])
wdh <- .3          # The width of the boxplots
scaler <- 20       # the ratio of the Fst and the dxy axis
clr_sec <- 'gray'  # the color of the secondary axis (dxy)

# start script -------------------

# import Fst
fst_files <- dir(fst_dir, pattern = '.50k.windowed.weir.fst.gz')

fst_data <- str_c(fst_dir,fst_files) %>%
  purrr::map(summarize_fst) %>%
  bind_rows()

# import dxy
dxy_files <- dir(dxy_dir)

dxy_data <-  str_c(dxy_dir,dxy_files) %>%
  purrr::map(summarize_dxy) %>%
  bind_rows()

fst_order <- fst_data %>%
  select(run, `mean_weighted-fst`) %>%
  mutate(run = fct_reorder(run,`mean_weighted-fst`))

data <- left_join(fst_data, dxy_data) %>%
  select(c(8,1:7,9:15)) %>%
  gather(key = 'stat', value = 'val',2:15) %>%
  separate(stat, into = c('sumstat','popstat'),sep = '_') %>%
  mutate(val_scaled = ifelse(popstat == 'dxy', val * scaler , val)) %>%
  unite(temp, val, val_scaled) %>%
  spread(.,key = 'sumstat',value = 'temp') %>%
  separate(mean, into = c('mean','mean_scaled'),sep = '_', convert = TRUE) %>%
  separate(median, into = c('median','median_scaled'),sep = '_', convert = TRUE) %>%
  separate(sd, into = c('sd','sd_scaled'),sep = '_', convert = TRUE) %>%
  separate(lower, into = c('lower','lower_scaled'),sep = '_', convert = TRUE) %>%
  separate(upper, into = c('upper','upper_scaled'),sep = '_', convert = TRUE) %>%
  separate(lowpoint, into = c('lowpoint','lowpoint_scaled'),sep = '_', convert = TRUE) %>%
  separate(highpoint, into = c('highpoint','highpoint_scaled'),sep = '_', convert = TRUE) %>%
  mutate(loc = str_sub(run,4,6),
         run = factor(run, levels = levels(fst_order$run)),
         x = as.numeric(run) ,
         x_dodge = ifelse(popstat == 'dxy',x + .25,x - .25),
         x_start_dodge = x_dodge - wdh/2,
         x_end_dodge = x_dodge + wdh/2,
         popstat_loc = str_c(popstat,'[',loc,']'))

run_ord <- tibble(run = levels(data$run),
                  run_ord = 1:length(levels(data$run)))

networx <- tibble( loc = c('bel','hon', 'pan'),
                   n = c(5,6,3),
                   label = list(str_c(c('ind','may','nig','pue','uni'),'bel'),
                                str_c(c('abe','gum','nig','pue','ran','uni'),'hon'),
                                str_c(c('nig','pue','uni'),'pan')),
                   weight = c(1,1.45,1)) %>%
  purrr::pmap(network_layout) %>%
  bind_rows()

plot_list <- networx %>%
  purrr::pmap(plot_network, node_lab_shift = .2)

p1 <- cowplot::plot_grid(
  grid::textGrob('Belize'),
  grid::textGrob('Honduras'),
  grid::textGrob('Panama'),
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  ncol = 3, rel_heights = c(.1,1))

p2 <- data %>%
  ggplot(aes(color = popstat_loc)) +
  geom_segment(aes(x = x_dodge, xend = x_dodge,
                   y = lowpoint_scaled,yend = highpoint_scaled))+
  geom_rect(aes(xmin = x_start_dodge, xmax = x_end_dodge,
                ymin = lower_scaled, ymax = upper_scaled),
             fill='white')+
  geom_segment(aes(x = x_start_dodge,
                   xend = x_end_dodge,
                   y = median_scaled,
                   yend = median_scaled),lwd = .9)+
  geom_point(aes(x = x_dodge, y = mean_scaled),shape = 21, size = .7, fill = 'white')+
  scale_x_continuous(breaks = 1:28) +
  scale_y_continuous(breaks = c(0,.05,.1,.15),
                     name = expression(italic(F[ST])),
                     sec.axis = sec_axis(~ . /scaler,
                                         name = expression(italic(d[XY])),
                                         breaks = c(0,.005,.01)))+
  scale_color_manual(values = c(make_faint_clr('bel'),
                                make_faint_clr('hon'),
                                make_faint_clr('pan'))[c(1,3,5,2,4,6)])+
  coord_cartesian(xlim = c(-1,29),
                  expand = c(0,0))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y.right = element_text(color = clr_sec),
        axis.title.y.right = element_text(color = clr_sec))

p_done <- cowplot::plot_grid(p1,p2,
                             ncol = 1,
                             rel_heights = c(.7,1),
                             labels = letters[1:2] %>% project_case())

hypo_save(p_done, filename = 'figures/F1.pdf',
          width = 9, height = 7,
          comment = plot_comment)