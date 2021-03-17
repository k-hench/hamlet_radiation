#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F1.R \
#    2_analysis/dxy/50k/ 2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/dxy/50k/', '2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt')
# script_name <- "R/fig/plot_F1.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)

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
run_ord <- tibble(run = levels(data$run),
                  run_ord = 1:length(levels(data$run)))

# underlying structure for the network plots
networx <- tibble( loc = c('bel','hon', 'pan'),
                   n = c(5,6,3),
                   label = list(str_c(c('ind','may','nig','pue','uni'),'bel'),
                                str_c(c('abe','gum','nig','pue','ran','uni'),'hon'),
                                str_c(c('nig','pue','uni'),'pan')),
                   weight = c(1,1.45,1)) %>%
  purrr::pmap(network_layout) %>%
  bind_rows() %>%
  mutate(edges = map(edges, function(x){x %>% left_join(data %>% filter(popstat == "weighted-fst") %>% select(run, median, mean)) }))


# plot the individual networks by location
plot_network <- function(loc, nodes, edges, asp = 0.8, sep = 0, node_lab_shift = 0){
  loc_edge <- c(bel = .68, hon = .66, pan = .82) -.03
  clr_prep <- (scales::colour_ramp(c("black",
                                     clr_loc[loc])))(c(0.4, 1))
  clrs <- colorRampPalette(clr_prep)(max(edges$idx))
  p <- nodes %>% ggplot(aes(x, y)) + 
    coord_fixed(ratio = asp) + 
    # ggfx:::with_blur.Layer(geom_arc(inherit.aes = FALSE,
    #          aes(x0 = x, 
    #              y0 = y, r = .3,
    #              start = -loc_edge[[loc]] * pi + seq(0,2*pi,length.out = max(idx+1))[idx], 
    #              end = loc_edge[[loc]] * pi + seq(0,2*pi,length.out = max(idx+1))[idx]),
    #          color = "black",#rgb(.8,.8,.8), 
    #          # arrow = arrow(type = "closed",length = unit(0,"pt")),
    #           linetype = 1, size = .6), 
    #          sigma = 1.15) +
    # geom_arc(inherit.aes = FALSE,
    #          aes(x0 = x, 
    #              y0 = y, r = .3,
    #              start = -loc_edge[[loc]] * pi + seq(0,2*pi,length.out = max(idx+1))[idx], 
    #              end = loc_edge[[loc]] * pi + seq(0,2*pi,length.out = max(idx+1))[idx],
    #              color = str_sub(label, -6,-4)), 
    #          linetype = 1, size = 1) +
    geom_segment(data = edges,
                 aes(xend = xend, yend = yend, size = median), 
                 color = clr_loc[loc]) +#, size = plot_lwd)
    scale_size(limits = c(0, 1), range = c(.1, 4))+
    scale_color_manual(values = clr)
  for (k in nodes$idx) {
    p <- p + plot_fish(short = str_sub(nodes$label[k], 1,  3),
                       x = nodes$x[k], y = nodes$y[k])
  }
  p + geom_label(data = edges, aes(x = xmid_shift + sign(xmid_shift) * sep,
                                   y = ymid_shift + sign(ymid_shift) * sep * asp,
                                   label = run_ord),
                 color = clr_loc[loc],
                 label.padding = unit(1, "pt"),
                 label.size = 0,
                 size = plot_text_size/ggplot2::.pt) +
    scale_fill_manual(values = clr_loc,
                      guide = FALSE) +
    scale_x_continuous(limits = c(-1.3, 1.3),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    theme_void()
}

plot_list <- networx %>%
    purrr::pmap(plot_network, node_lab_shift = .2)

pca_plot <- function(loc){
  evs <- str_c("2_analysis/pca/", loc ,".exp_var.txt.gz") %>% 
    read_tsv()
  str_c("2_analysis/pca/", loc ,".scores.txt.gz") %>% 
    read_tsv() %>% 
    mutate(spec = str_sub(id, -6,-4)) %>%
    ggplot(aes(x = EV01, y = EV02, fill = spec))+
    geom_point(shape = 21) +
    labs(x = str_c("EV01 (", sprintf("%.2f",evs$exp_var[[1]]), " %)"),
         y = str_c("EV02 (", sprintf("%.2f",evs$exp_var[[2]]), " %)"))+
    scale_fill_manual(values = clr)+
    theme_minimal()+
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          panel.background = element_rect(color = "transparent", fill = "#eeeeee"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = plot_text_size))
}

pcas <- c("bel", "hon", "pan") %>% map(pca_plot)

# assemble panel a
p1 <- cowplot::plot_grid(
  grid::textGrob('Belize', gp = gpar(fontsize = plot_text_size)),
  grid::textGrob('Honduras', gp = gpar(fontsize = plot_text_size)),
  grid::textGrob('Panama', gp = gpar(fontsize = plot_text_size)),
  plot_list[[1]] + theme(legend.position = "none"), plot_list[[2]] + theme(legend.position = "none"), plot_list[[3]] + theme(legend.position = "none"),
  pcas[[1]], pcas[[2]], pcas[[3]],
  # NULL, 
  # get_legend(plot_list[[1]] + guides(size = guide_legend(title = expression(Median~italic(F[ST])),title.position = "top",
  #                                                        override.aes = list(color = "#aaaaaa"))) + theme(legend.position = "bottom")),
  # NULL,
  ncol = 3, rel_heights = c(.1,1,1),
  labels = c(letters[1], "", "", "", "","", letters[2:4]) %>%
    project_case(),
  label_size = plot_text_size)

# assemble panel b
p2 <- data %>%
  filter(popstat == "weighted-fst") %>%
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
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
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
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y.right = element_text(color = clr_sec),
        axis.title.y.right = element_text(color = clr_sec))

# merge panel b & c
p23 <- cowplot::plot_grid(p2,p3,
                          ncol = 2,
                          labels = letters[5:6] %>%
                            project_case(),
                          label_size = plot_text_size)

# merge all panels
p_done <- cowplot::plot_grid(p1, p23,
                             ncol = 1,
                             rel_heights = c(1.6,1))
# export figure 1
hypo_save(p_done, filename = 'figures/F1_expand.pdf',
          width = f_width, height =  f_width * 1,
          comment = plot_comment)

# compile fst and dxy table (table 1) for the manuscript
table_all <- dxy_data %>%
  select(run, mean_dxy) %>%
  left_join( vroom::vroom(fst_globals, delim = '\t',
                          col_names = c('loc','run','mean','weighted_fst')) %>%
               mutate(run = str_c(loc,'-',run) %>%
                        reformat_run_name())  %>%
               select(run, weighted_fst)) %>%
  pivot_longer(names_to = 'stat',2:3) %>%
  separate(run, into = c('pop1', 'pop2'), sep = '-') %>%
  mutate(prep1 = ifelse(stat == "weighted_fst", pop2,pop1),
         prep2 = ifelse(stat == "weighted_fst", pop1,pop2),
         pop1 = factor(prep1, levels = pop_levels),
         pop2 = factor(prep2, levels = pop_levels),
         value = sprintf('%7.5f', value) ) %>%
  select(pop1,pop2,value) %>%
  arrange(pop2,pop1) %>%
  mutate(pop2 = as.character(pop2) %>%
           str_replace(pattern = '([a-z]{3})([a-z]{3})',
                       replacement = '\\1|\\2'),
         pop1 = as.character(pop1) %>%
           str_replace(pattern = '([a-z]{3})([a-z]{3})',
                       replacement = '\\1|\\2')) %>%
  pivot_wider(values_from = value,
              names_from = pop2) %>%
  rename( Population = 'pop1') %>%
  mutate(srt1 = str_sub(Population,-3, -1),
         srt2 = str_sub(Population,1, 3))  %>%
  arrange(srt1,srt2) %>%
  select(-srt1,-srt2)

# replace "NA" by dashes "-"
table_all[is.na(table_all)] <- '-'

# export sub-tables 1 a - c
table_all[1:5,c(1,2:6)] %>% export_2_latex(name = 'tables/suppl_tab3a.tex')
table_all[6:11,c(1,7:12)] %>% export_2_latex(name = 'tables/suppl_tab3b.tex')
table_all[12:14, c(1,13:15)] %>% export_2_latex(name = 'tables/suppl_tab3c.tex')
