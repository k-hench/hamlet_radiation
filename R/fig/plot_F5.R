#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F5.R \
#   2_analysis/twisst/weights/ ressources/plugin/trees/ \
#   https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
#   2_analysis/summaries/fst_outliers_998.tsv 2_analysis/dxy/50k/ \
#   2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt \
#   2_analysis/GxP/50000/ 200 5 2_analysis/fst/poptree/summary/
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
#           "2_analysis/fst/poptree/summary/")
# script_name <- "R/fig/plot_F5.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(furrr)
library(ggraph)
library(tidygraph)
library(ggtext)
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
source(twisst_functions, local = TRUE)

plan(multiprocess)
window_buffer <- 2.5*10^5
#-------------------
library(ape)
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

# compose base figure
p_single <- outlier_table %>%
  filter(outlier_id %in% outlier_pick) %>%
  left_join(region_label_tibbles) %>%
  mutate(outlier_nr = row_number(),
         text = ifelse(outlier_nr == 1,TRUE,FALSE),
         trait = c('Snout', 'Bars', 'Peduncle')) %>%
  pmap(plot_curtain, cool_genes = cool_genes) %>%
  c(., poptree_plot_list) %>%
  cowplot::plot_grid(plotlist = ., nrow = 2,
                     rel_heights = c(1, .17),
                     labels = letters[1:length(outlier_pick)] %>% project_case())

# compile legend
# dummy plot for fst legend
p_dummy_fst <- outlier_table %>% filter(row_number() == 1) %>% purrr::pmap(plot_panel_fst) %>% .[[1]]
# dummy plot for  gxp legend
p_dummy_gxp <- outlier_table %>% filter(row_number() == 1) %>% purrr::pmap(plot_panel_gxp, trait = 'Bars') %>% .[[1]]
# fst legend
p_leg_fst <- (p_dummy_fst+theme(legend.position = 'bottom')) %>% get_legend()
# gxp legend
p_leg_gxp <- (p_dummy_gxp+theme(legend.position = 'bottom')) %>% get_legend()
# create poptree legend
p_leg_poptree <- (poptree_plot_list[[1]] + theme(legend.position = "bottom")) %>% get_legend()
# create sub-legend 1
p_leg1 <- cowplot::plot_grid(p_leg_fst,
                             p_leg_gxp,
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
                            nrow = 1, rel_widths = c(.6, 1))
# complete legend
p_leg <- cowplot::plot_grid(
  p_leg_poptree,
  p_leg_3,
  ncol = 1,
  rel_heights = c(.2, 1))

# finalize figure
p_done <- cowplot::plot_grid(p_single, p_leg,
                             ncol = 1,
                             rel_heights = c(1, .17))

# export figure
scl <- .9
hypo_save(plot = p_done, filename = 'figures/F5.pdf',
          width = 14*scl, height = 13*scl,
          comment = script_name,
          device = cairo_pdf)
