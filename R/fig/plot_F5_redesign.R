#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F5_redesign.R \
#   2_analysis/twisst/weights/ ressources/plugin/trees/ \
#   https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
#   2_analysis/summaries/fst_outliers_998.tsv 2_analysis/dxy/50k/ \
#   2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt \
#   2_analysis/GxP/50000/ 200 5 2_analysis/fst/poptree/summary/ \
#   ~/work/puebla_lab/stash/pomo/
# ===============================================================
# This script produces Figure 5 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/twisst/weights/', 'ressources/plugin/trees/',
#           'https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R',
#           '2_analysis/summaries/fst_outliers_998.tsv',
#           '2_analysis/dxy/50k/', '2_analysis/fst/50k/',
#           '2_analysis/summaries/fst_globals.txt',
#           '2_analysis/GxP/50000/', 200, 5,
#           "2_analysis/fst/poptree/summary/",
#           "~/work/puebla_lab/stash/pomo/")
# script_name <- "R/fig/plot_F5.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(furrr)
library(ggraph)
library(tidygraph)
library(ggtext)
library(ape)
library(ggtree)
library(patchwork)
library(ggplotify)
library(phangorn)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
w_path <- as.character(args[1])
d_path <- as.character(args[2])
twisst_functions <- as.character(args[3])
out_table <- as.character(args[4])
dxy_dir <- as.character(args[5])
fst_dir <- as.character(args[6])
fst_globals <- as.character(args[7])
gxp_dir <- as.character(args[8])
twisst_size <- as.numeric(args[9])
resolution <- as.numeric(args[10])
fst_summary_path <- as.character(args[11])

pomo_path <- as.character(args[12])
pomo_trees <- dir(pomo_path, pattern = "155_pop")

source(twisst_functions, local = TRUE)

plan(multiprocess)
window_buffer <- 2.5*10^5
#-------------------
library(igraph)
# actual script =========================================================
# locate dxy data files
dxy_files <- dir(dxy_dir, pattern = str_c('dxy.*[a-z]{3}.*.', resolution ,'0kb-', resolution ,'kb.tsv.gz'))

# import dxy data
dxy_data <- tibble(file = str_c(dxy_dir, dxy_files)) %>%
  purrr::pmap_dfr(get_dxy, kb = str_c(resolution ,'0kb'))

# set traits of interest for GxP
gxp_traits <- c('Bars', 'Snout', 'Peduncle')

# import GxP data
gxp_data <- str_c(gxp_dir,gxp_traits,'.lm.', resolution ,'0k.', resolution ,'k.txt.gz') %>%
  future_map_dfr(get_gxp_long, kb = 50)

# set topology weighting color scheme
twisst_clr <- c(Blue = "#0140E5", Bars = "#E32210", Butter = "#E4E42E")
# set GxP color scheme
gxp_clr <- c(Bars = "#79009f", Snout = "#E48A00", Peduncle = "#5B9E2D") %>%
  darken(factor = .95) %>%
  set_names(., nm = gxp_traits)

# copute genome wide average dxy
dxy_globals <- dxy_data %>%
  filter(BIN_START %% ( resolution * 10000 ) == 1 ) %>%
  group_by( run ) %>%
  summarise(mean_global_dxy = sum(dxy*N_SITES)/sum(N_SITES)) %>%
  mutate(run = fct_reorder(run,mean_global_dxy))

# import genome wide average fst
# and order population pairs by genomewide average fst
fst_globals <- vroom::vroom(fst_globals,delim = '\t',
                        col_names = c('loc','run_prep','mean_fst','weighted_fst')) %>%
  separate(run_prep,into = c('pop1','pop2'),sep = '-') %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted_fst))

# locate fst data files
fst_files <- dir(fst_dir ,pattern = str_c('.', resolution ,'0k.windowed.weir.fst.gz'))
# import fst data
fst_data <- str_c(fst_dir, fst_files) %>%
  future_map_dfr(get_fst, kb = str_c(resolution ,'0k')) %>%
  left_join(dxy_globals) %>%
  left_join(fst_globals) %>%
  mutate(run = refactor(., fst_globals),
         BIN_MID = (BIN_START+BIN_END)/2)

# order dxy population pairs by genomewide average fst
dxy_data <- dxy_data %>%
  left_join(dxy_globals) %>%
  left_join(fst_globals) %>%
  mutate(run = refactor(dxy_data, fst_globals),
         window = 'bold(italic(d[xy]))')

# compute delta dxy
data_dxy_summary <- dxy_data %>%
  group_by(GPOS) %>%
  summarise(scaffold = CHROM[1],
            start = BIN_START[1],
            end = BIN_END[1],
            mid = BIN_MID[1],
            min_dxy = min(dxy),
            max_dxy = max(dxy),
            mean_dxy = mean(dxy),
            median_dxy = median(dxy),
            sd_dxy = sd(dxy),
            delta_dxy = max(dxy)-min(dxy))

# twisst part ------------------
# load fst outlier regions
outlier_table <- vroom::vroom(out_table, delim = '\t') %>%
  setNames(., nm = c("outlier_id","lg", "start", "end", "gstart","gend","gpos"))

# set outlier regions of interest
outlier_pick = c('LG04_1', 'LG12_3', 'LG12_4')

# select genes to label
cool_genes <-  c('arl3','kif16b','cdx1','hmcn2',
                 'sox10','smarca4',
                 'rorb',
                 'alox12b','egr1',
                 'ube4b','casz1',
                 'hoxc8a','hoxc9','hoxc10a',
                 'hoxc13a','rarga','rarg',
                 'snai1','fam83d','mafb','sws2abeta','sws2aalpha','sws2b','lws','grm8')

# load twisst data
data_tables <- list(bel = prep_data(loc = 'bel'),
                    hon = prep_data(loc = 'hon'))

# set species sampled in belize
pops_bel <- c('ind', 'may', 'nig', 'pue', 'uni')

# set outlier region label
region_label_tibbles <- tibble(outlier_id = outlier_pick,
                            label = c('A','B','C'))

# load and recolor trait icons
trait_grob <- tibble(svg = hypoimg::hypo_trait_img$grob_circle[hypoimg::hypo_trait_img$trait %in% gxp_traits],
                     layer = c(4,3,7),
                     color = gxp_clr[gxp_traits %>% sort()])%>%
  pmap(.f = hypo_recolor_svg) %>%
  set_names(nm = gxp_traits %>% sort())

# recolor second bars-layer
trait_grob[["Bars"]] <- trait_grob[["Bars"]] %>% hypo_recolor_svg(layer = 7,color = gxp_clr[["Bars"]])

# prepare population tree trait_panels ----------------------
# load outler region wide average fst data
tree_list <- outlier_pick %>%
  purrr::map_dfr(get_fst_summary_data)

# create neighbour-joining trees and plot
poptree_plot_list <-  tree_list %>%
  purrr::pmap(plot_fst_poptree)

import_tree <- function(file){
  tag = file %>%
    stringr::str_remove(".*/") %>%
    stringr::str_remove("_pop.cf.treefile")
  tibble::tibble(outlier = stringr::str_sub(tag, 1, 6) %>% 
           stringr::str_to_upper() %>% 
           stringr::str_replace("\\.","_"),
         type = stringr::str_sub(tag, 8, -1) ,
         tree = list(read.tree(file)))
}

root_hamlets <- function(tree, type){
  outgr <- switch (type,
                   "155" = "floflo",
                   "hyS" = c("tabhon", "torpan")
  )
  ape::root(phy = tree, outgroup = outgr)  
}

conditional_rotate <- function(tree, rotate = TRUE, root_node){
  if(rotate){tree %>% ggtree::rotate(root_node)} else{tree}
}

conditional_highlight <- function(tree,
                                  highl = TRUE, higl_node, 
                                  support_guide = FALSE ){
  if(highl){
    p <- tree %>% 
      ggtree::groupClade(.node = higl_node, group_name =  "clade") %>% 
      ggtree(layout = "circular", aes(color = clade), size = plot_lwd) 
  } else {
    p <- tree %>% 
      ggtree(layout = "circular", size = plot_lwd)
  }
  p <- p +
    geom_nodepoint(aes(fill =  support_class, size = support_class),
                   shape = 21) +
    scale_fill_manual(values = c(`(0,0.5]` = "transparent",
                                 `(0.5,0.7]` = "white",
                                 `(0.7,0.9]` = "gray",
                                 `(0.9,1]` = "black"),
                      drop = FALSE,
                      guide = FALSE) +
    scale_size_manual(values = c(`(0,0.5]` = 0,
                                 `(0.5,0.7]` = .9,
                                 `(0.7,0.9]` = .9,
                                 `(0.9,1]` = .9),
                      na.value = 0,
                      drop = FALSE,
                      guide = FALSE) 
  
  if(support_guide){p <- p +   guides(fill = guide_legend(title = "Support", row = 1,keywidth = unit(5, "pt")),
                                      size = guide_legend(title = "Support", row = 1,keywidth = unit(5, "pt"))) }
  p
}

plot_tree <- function(tree, angle_in = 0,
                      color = "red", 
                      higl_clade = TRUE, higl_node = NA,
                      root_node = NA, rotate_root = FALSE, 
                      xlim = c(-2,4.5), 
                      support_guide = FALSE,
                      open_angle = 180){
  p <- (tree %>% 
     ggtree(layout = "circular") %>% 
     .$data %>% 
     mutate(support = as.numeric(label) /100,
            support_class = cut(support, c(0,.5,.7,.9,1)) %>% 
              as.character() %>%
              factor(levels = c("(0,0.5]", "(0.5,0.7]", "(0.7,0.9]", "(0.9,1]"))
     ) %>% 
     conditional_highlight(highl = higl_clade,higl_node = higl_node, support_guide = support_guide) %>% 
     conditional_rotate(rotate = rotate_root, root_node = root_node) %>% 
     open_tree(angle =  open_angle) %>% 
     rotate_tree(angle = angle_in)) +
    geom_tiplab(offset = diff(xlim) * .07, size = plot_text_size / ggplot2::.pt *.7#, 
                # aes(label = str_c(str_sub(label, -6,-6),str_c(str_sub(label, -3,-3))) %>% 
                #       str_to_upper())
                ) +
    # geom_tippoint(size = .4) +
    scale_color_manual(values  = c(`0` = "black", `1` = color),
                       guide = FALSE) +
    scale_x_continuous(limits = xlim) +
    theme_void()
  
  ggplot() +
    coord_equal(xlim = c(-.1, .95), 
                ylim = c(-.01, .52),
                expand = 0) +
    # annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Snout"][[1]],
    #                   xmin = .4, xmax = .6,
    #                   ymin = .4, ymax = .6)+
    annotation_custom(grob = ggplotGrob(p),
    ymin = 0 , ymax = 1,
    xmin = 0, xmax = 1) +
    theme_void()
}

pomo_data <- str_c(pomo_path, pomo_trees) %>% 
  purrr::map_dfr(import_tree) %>% 
  mutate(rooted_tree = map2(.x = tree,.y = type, .f = root_hamlets))

p_pomo1 <- plot_tree(midpoint(pomo_data$tree[[1]]), 
          higl_node = 21,
          angle_in = 168,
          color = twisst_clr[["Blue"]], 
          xlim = c(-.015, .1),
          open_angle = 168)

p_pomo2 <- plot_tree(midpoint(pomo_data$tree[[2]]), 
          higl_node = 27,
          color = twisst_clr[["Bars"]], 
          xlim = c(-.015,.065),
          angle_in = 168,
          open_angle = 168)

p_pomo3 <- plot_tree(midpoint(pomo_data$tree[[3]]), 
          higl_node = 28,
          color = twisst_clr[["Butter"]], 
          xlim = c(-.01, .053),
          angle_in = 168,
          open_angle = 168)

plot_list <- list(p_pomo1, p_pomo2, p_pomo3)

# compose base figure
p_single <- outlier_table %>%
  filter(outlier_id %in% outlier_pick) %>%
  left_join(region_label_tibbles) %>%
  mutate(outlier_nr = row_number(),
         text = ifelse(outlier_nr == 1,TRUE,FALSE),
         trait = c('Snout', 'Bars', 'Peduncle')) %>%
  pmap(plot_curtain, cool_genes = cool_genes) %>%
   c(., plot_list) %>%
  cowplot::plot_grid(plotlist = ., nrow = 2,
                     rel_heights = c(1, .18),
                     labels = letters[1:length(outlier_pick)] %>% project_case(),
                     label_size = plot_text_size)

# compile legend
# dummy plot for fst legend
p_dummy_fst <- outlier_table %>% filter(row_number() == 1) %>% 
  purrr::pmap(plot_panel_fst) %>% 
  .[[1]]+ guides(color = guide_colorbar(barheight = unit(3, "pt"),
                                        barwidth = unit(100, "pt"),
                                        label.position = "top",
                                        ticks.colour = "black"))
# dummy plot for  gxp legend
p_dummy_gxp <- outlier_table %>% filter(row_number() == 1) %>% purrr::pmap(plot_panel_gxp, trait = 'Bars') %>% .[[1]]
# fst legend
p_leg_fst <- (p_dummy_fst+theme(legend.position = 'bottom')) %>% get_legend()
# gxp legend
p_leg_gxp <- (p_dummy_gxp+theme(legend.position = 'bottom')) %>% get_legend()
# create poptree legend
# p_leg_poptree <- (poptree_plot_list[[1]] + 
#                     theme(text = element_text(size = plot_text_size),
#                           legend.position = "bottom")) %>% get_legend()
# create sub-legend 1

p_leg_pomo <- ((midpoint(pomo_data$tree[[1]]) %>% 
    ggtree(layout = "circular") %>% 
    .$data %>% 
    mutate(support = as.numeric(label) /100,
           support_class = cut(support, c(0,.5,.7,.9,1)) %>% 
             as.character() %>%
             factor(levels = c("(0,0.5]", "(0.5,0.7]", "(0.7,0.9]", "(0.9,1]")))) %>% 
conditional_highlight(tree = ., 
          higl_node = 21,
          highl = FALSE,
          support_guide = TRUE) +
  theme(text = element_text(size = plot_text_size),
        legend.position = "bottom") ) %>%
  get_legend()

p_leg1 <- cowplot::plot_grid(p_leg_fst,
                             p_leg_gxp,
                             p_leg_pomo,
                             ncol = 1)

# create sub-legend 2 (phylogeny schematics)
p_leg2 <- tibble(spec1 = c('indigo', 'indigo','unicolor'),
                 spec2 = c('maya', 'puella',NA),
                 color = twisst_clr %>% unname() %>% darken(.,factor = .8),
                 mode = c(rep('pair',2),'isolation')) %>%
  future_pmap(plot_leg) %>%
  cowplot::plot_grid(plotlist = .,
                     nrow = 1)

# create sub-legend 3
p_leg_3 <- cowplot::plot_grid(p_leg1,
                            p_leg2,
                            nrow = 1, 
                            rel_widths = c(.6, 1))
# complete legend
# p_leg <- cowplot::plot_grid(
#   p_leg_poptree,
#   p_leg_3,
#   ncol = 1,
#   rel_heights = c(.2, 1))

# finalize figure
p_done <- cowplot::plot_grid(p_single, p_leg_3,
                             ncol = 1,
                             rel_heights = c(1, .17))

# export figure
hypo_save(plot = p_done, filename = 'figures/F5_redesign.pdf',
          width = f_width, 
          height = f_width * .93,
          comment = script_name,
          device = cairo_pdf)
  
