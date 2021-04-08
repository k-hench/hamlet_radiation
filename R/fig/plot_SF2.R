#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF2.R \
#    2_analysis/dxy/50k/ 2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Suppl. Figure 2 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/dxy/50k/', '2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt')
# script_name <- "R/fig/plot_SF2.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(patchwork)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
dxy_dir <- as.character(args[1])
fst_dir <- as.character(args[2])
fst_globals <- as.character(args[3])
wdh <- .3          # The width of the boxplots
scaler <- 20       # the ratio of the Fst and the dxy axis
clr_sec <- 'gray'  # the color of the secondary axis (dxy)

# start script -------------------

# import Fst
fst_files <- dir(fst_dir, pattern = '.50k.windowed.weir.fst.gz')

fst_data <- str_c(fst_dir,fst_files) %>%
  purrr::map(summarize_fst) %>%
  bind_rows()

# lookup dxy files
dxy_files <- dir(dxy_dir)

# import dxy
dxy_data <-  str_c(dxy_dir,dxy_files) %>%
  purrr::map(summarize_dxy) %>%
  bind_rows()

# determine fst ranking
fst_order <- fst_data %>%
  select(run, `mean_weighted-fst`) %>%
  mutate(run = fct_reorder(run, `mean_weighted-fst`))

# merge fst and dxy cc_data
# (large parts of this code are now unnecessary after the separation of dxy and
#  fst plots into separate panels b & c)
data <- left_join(fst_data, dxy_data) %>%
  select(c(8,1:7,9:15)) %>%
  # reformat table to enable parallel plotting (with secondary axis)
  gather(key = 'stat', value = 'val', 2:15) %>%
  # sumstat contains the values needed to plot the boxplots (quartiles, etc)
  separate(stat, into = c('sumstat', 'popstat'), sep = '_') %>%
  # duplicate dxy values scaled to fst range
  mutate(val_scaled = ifelse(popstat == 'dxy', val * scaler , val)) %>%
  unite(temp, val, val_scaled) %>%
  # separate th eoriginal values from the scales ons (scaled = secondary axis)
  spread(.,key = 'sumstat',value = 'temp') %>%
  separate(mean, into = c('mean','mean_scaled'),sep = '_', convert = TRUE) %>%
  separate(median, into = c('median','median_scaled'), sep = '_', convert = TRUE) %>%
  separate(sd, into = c('sd','sd_scaled'),sep = '_', convert = TRUE) %>%
  separate(lower, into = c('lower','lower_scaled'), sep = '_', convert = TRUE) %>%
  separate(upper, into = c('upper','upper_scaled'), sep = '_', convert = TRUE) %>%
  separate(lowpoint, into = c('lowpoint','lowpoint_scaled'), sep = '_', convert = TRUE) %>%
  separate(highpoint, into = c('highpoint','highpoint_scaled'), sep = '_', convert = TRUE) %>%
  # include "dodge"-positions for side-by-side plotting (secondary axis)
  mutate(loc = str_sub(run,4,6),
         run = factor(run, levels = levels(fst_order$run)),
         x = as.numeric(run) ,
         x_dodge = ifelse(popstat == 'dxy', x + .25, x - .25),
         x_start_dodge = x_dodge - wdh/2,
         x_end_dodge = x_dodge + wdh/2,
         popstat_loc = str_c(popstat,'[',loc,']'))

# sort run by average genome wide Fst


fst_data_gather <- data %>%
  filter(popstat == "weighted-fst")  %>% 
  gather(key = 'stat', value = 'val', -run) %>%
  # sumstat contains the values needed to plot the boxplots (quartiles, etc)
  separate(stat, into = c('sumstat', 'popstat'), sep = '_') %>%
  # duplicate dxy values scaled to fst range
  mutate(val_scaled = ifelse(popstat == 'dxy', val * scaler , val)) %>%
  unite(temp, val, val_scaled) %>%
  # separate th eoriginal values from the scales ons (scaled = secondary axis)
  spread(.,key = 'sumstat',value = 'temp') %>%
  separate(mean, into = c('mean','mean_scaled'),sep = '_', convert = TRUE) %>%
  separate(median, into = c('median','median_scaled'), sep = '_', convert = TRUE) %>%
  separate(sd, into = c('sd','sd_scaled'),sep = '_', convert = TRUE) %>%
  separate(lower, into = c('lower','lower_scaled'), sep = '_', convert = TRUE) %>%
  separate(upper, into = c('upper','upper_scaled'), sep = '_', convert = TRUE) %>%
  separate(lowpoint, into = c('lowpoint','lowpoint_scaled'), sep = '_', convert = TRUE) %>%
  separate(highpoint, into = c('highpoint','highpoint_scaled'), sep = '_', convert = TRUE) %>%
  # include "dodge"-positions for side-by-side plotting (secondary axis)
  mutate(loc = str_sub(run,4,6),
         run = factor(run, levels = levels(fst_order$run)),
         x = as.numeric(run) ,
         x_dodge = ifelse(popstat == 'dxy', x + .25, x - .25),
         x_start_dodge = x_dodge - wdh/2,
         x_end_dodge = x_dodge + wdh/2,
         popstat_loc = str_c(popstat,'[',loc,']'))

run_ord <- tibble(run = levels(data$run),
                  run_ord = 1:length(levels(data$run)))


networx <- tibble( loc = c('bel','hon', 'pan'),
                   n = c(5, 6, 3),
                   label = list(str_c(c('ind','may','nig','pue','uni'),'bel'),
                                str_c(c('abe','gum','nig','pue','ran','uni'),'hon'),
                                str_c(c('nig','pue','uni'),'pan')),
                   weight = c(1,1.45,1)) %>%
  purrr::pmap_dfr(network_layout) %>%
  mutate(edges = map(edges, function(x){x %>% left_join(fst_data_gather %>% filter(popstat == "weighted-fst") %>% select(run, median, mean)) }))


hypo_anno_l_lwd <- function (species, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, lwd = .15, line_color = "black") {
  stopifnot(length(species) == 1)
  stopifnot(is.character(species))
  stopifnot(species %in% hypo_img$spec)
  nr_species <- which(hypo_img$spec == species)
  annotation_custom(editGrob(grob = hypo_img$l[[nr_species]],
                             gPath = "GRID.picComplexPath.*", grep = TRUE,
                             gp = gpar( lwd = lwd, col = line_color
                             ), 
                             global = TRUE, strict = FALSE) ,
                    xmin = xmin, 
                    xmax = xmax, ymin = ymin, ymax = ymax)
}

plot_fish_lwd <- function (short, x = 0, y = 3, height = 5, width = 5, lwd = .15, line_color = "transparent") { 
  hypo_anno_l_lwd(sp_names[short], xmin = x - 0.5 * width, xmax = x + 
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 * 
                    height, lwd = lwd, line_color = line_color)
}


plot_network <- function(loc, nodes, edges, asp = 0.8, sep = 0, node_lab_shift = 0){
  loc_edge <- c(bel = .68, hon = .66, pan = .82) -.03
  clr_prep <- (scales::colour_ramp(c("black",
                                     clr_loc[loc])))(c(0.4, 1))
  clrs <- colorRampPalette(clr_prep)(max(edges$idx))
  p <- nodes %>% ggplot(aes(x, y)) + 
    coord_fixed(ratio = asp) + 
    geom_segment(data = edges,
                 aes(xend = xend, yend = yend),#, size = median), 
                 size = .1,
                 color = clr_loc[loc]) +#, size = plot_lwd)
    # scale_size(limits = c(0, 1), range = c(.1, 4))+
    scale_size(limits = c(0, 1), range = c(.1, 2))+
    scale_color_manual(values = clr)
  for (k in nodes$idx) {
    p <- p + plot_fish_lwd(short = str_sub(nodes$label[k], 1,  3),
                           x = nodes$x[k], y = nodes$y[k], height = .7, width = .7)
  }
  p + geom_label(data = edges, aes(x = xmid_shift + sign(xmid_shift) * sep,
                                   y = ymid_shift + sign(ymid_shift) * sep * asp,
                                   label = run_ord),
                 color = clr_loc[loc],
                 label.padding = unit(1, "pt"),
                 label.size = 0,
                 size = plot_text_size * .5 /ggplot2::.pt) +
    scale_fill_manual(values = clr_loc,
                      guide = FALSE) +
    scale_x_continuous(limits = c(-1.3, 1.3),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    theme_void()
}


plot_list <- networx %>%
  purrr::pmap(plot_network, node_lab_shift = .2)

p_net <- cowplot::plot_grid(
  # grid::textGrob('Belize', gp = gpar(fontsize = plot_text_size,)),
  # grid::textGrob('Honduras', gp = gpar(fontsize = plot_text_size)),
  # grid::textGrob('Panama', gp = gpar(fontsize = plot_text_size)),
  plot_list[[1]] + theme(legend.position = "none"), plot_list[[2]] + theme(legend.position = "none"), plot_list[[3]] + theme(legend.position = "none"),
  ncol = 3#, 
  #rel_heights = c(.1,1)
) %>% cowplot::as_grob()


# assemble panel b
p2 <- data %>%
  filter(popstat == "weighted-fst") %>%
  ggplot(aes(color = loc)) + 
  annotation_custom(p_net, ymin = .15, xmax = 25) +
  geom_segment(aes(x = x, xend = x,
                   y = lowpoint, yend = highpoint),
               lwd = plot_lwd)+
  geom_rect(aes(xmin = x - wdh, xmax = x + wdh,
                ymin = lower, ymax = upper),
            fill = 'white',
            size = plot_lwd)+
  geom_segment(aes(x = x - wdh,
                   xend = x + wdh,
                   y = median,
                   yend = median),
               lwd = plot_lwd)+
  geom_point(aes(x = x, y = mean),
             shape = 21, size = .7, fill = 'white')+
  scale_x_continuous(breaks = 1:28) +
  scale_y_continuous(#breaks = c(0,.05,.1,.15),
    name = expression(italic(F[ST])))+
  scale_color_manual(values = c(make_faint_clr('bel'),
                                make_faint_clr('hon'),
                                make_faint_clr('pan'))[c(2,4,6)])+
  coord_cartesian(xlim = c(0,29),
                  expand = c(0,0))+
  theme_minimal()+
  theme(text = element_text(size = plot_text_size),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        axis.text.y.right = element_text(color = clr_sec),
        axis.title.y.right = element_text(color = clr_sec))

# assemble panel c
p3 <- data %>%
  filter(popstat == "dxy") %>%
  ggplot(aes(color = loc)) +
  geom_segment(aes(x = x, xend = x,
                   y = lowpoint, yend = highpoint),
               lwd = plot_lwd)+
  geom_rect(aes(xmin = x - wdh, xmax = x + wdh,
                ymin = lower, ymax = upper),
            fill = 'white',
            size = plot_lwd)+
  geom_segment(aes(x = x - wdh,
                   xend = x + wdh,
                   y = median,
                   yend = median),
               lwd = plot_lwd)+
  geom_point(aes(x = x, y = mean),
             shape = 21, size = .7, fill = 'white')+
  scale_x_continuous(breaks = 1:28) +
  scale_y_continuous( expression(italic(d[XY])),
                      breaks = c(0,.0025,.005,.0075,.01),
                      limits = c(0,.01))+
  scale_color_manual(values = c(make_faint_clr('bel'),
                                make_faint_clr('hon'),
                                make_faint_clr('pan'))[c(2,4,6)])+
  coord_cartesian(xlim = c(0,29),
                  expand = c(0,0))+
  theme_minimal()+
  theme(text = element_text(size = plot_text_size),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        axis.text.y.right = element_text(color = clr_sec),
        axis.title.y.right = element_text(color = clr_sec))

# merge panel b & c
p_done <- p2 + p3 + plot_annotation(tag_levels = "a")

# export figure 1
hypo_save(p_done, filename = 'figures/SFY5.pdf',
          width = f_width, 
          height =  f_width * .4,
          comment = plot_comment)

hypo_save(plot_grid(p_net, p3, ncol = 1,
                    rel_heights = c(.4, 1),
                    labels = letters[1:2],
                    label_fontface = "plain",
                    label_size = plot_text_size), 
          filename = 'figures/SF2.pdf',
          width = f_width * .5, 
          height =  f_width * .5,
          comment = plot_comment)
