#!/usr/bin/env Rscript
#
# Context: this script depends on the input file 2_analysis/summaries/fst_permutation_summary.tsv
#          which is created by R/figs/plot_SFY.R
#
# run from terminal:
# Rscript --vanilla R/fig/plot_F1_redesign.R \
#    2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt 2_analysis/summaries/fst_permutation_summary.tsv
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt', '2_analysis/summaries/fst_permutation_summary.tsv')
# script_name <- "R/fig/plot_F1.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(patchwork)
library(ape)
library(ggraph)
library(tidygraph)
library(stringr)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
# config -----------------------
fst_dir <- as.character(args[1])
fst_globals <- as.character(args[2])
fst_permutation_file <- as.character(args[3])
wdh <- .3          # The width of the boxplots
scaler <- 20       # the ratio of the Fst and the dxy axis (legacy - not really needed anymore)
clr_sec <- 'gray'  # the color of the secondary axis (dxy)
# start script -------------------
# === fig 1. panel a: modified plugin from ressources/Rabosky_etal_2018/scripts/main figures/Figure3.R =============
library(BAMMtools)
library(geiger)
library(ggplotify)
source("R/bammtools_plot_mod.R")

basepath <- 'ressources/Rabosky_etal_2018/'
source(paste0(basepath, "scripts/supporting_fxns/PlottingFunctions.R"))

eventdata_vr <- paste0(basepath, "dataFiles/bamm_results/12k_tv1/event_data_thinned.csv")
eventdata_cr <- paste0(basepath, "dataFiles/bamm_results/12k_tc1/event_data_thinned.csv")
treefile <- paste0(basepath, "dataFiles/bamm_results/12k_tv1/bigfish_no_outgroup.tre")
fspdata <- paste0(basepath, "dataFiles/rate_lat_stats_by_sp_fixed0.5.csv")

anadromous <- FALSE

vx <- read.tree(treefile) 

# node rotations for plotting:
rset <- c(15447, 15708:15719)
for (i in 1:length(rset)){
  vx <- rotate(vx, node = rset[i])
}

vx <- read.tree(text = write.tree(vx))

spdata <- read.csv(fspdata, stringsAsFactors = F)
rownames(spdata) <- spdata$sp

# ---------------------------

latvals <- abs(spdata$lat_centroid)
names(latvals) <- spdata$sp

inboth <- intersect(spdata$sp, vx$tip.label)

# 1 species is in tree, but was dropped because is synonym to another species in tree according to fishbase (Gadus_ogac matches to Gadus_macrocephalus)
inboth <- intersect(inboth, spdata$sp[which(!is.na(spdata$tv.lambda))]) 

latvals <- latvals[inboth]

edvr <- getEventData(vx, eventdata_vr, burnin=0)

serranids <- c("Hypoplectrus_gemma", "Hypoplectrus_unicolor", "Hypoplectrus_gummigutta",
               "Hypoplectrus_chlorurus", "Hypoplectrus_aberrans", "Hypoplectrus_nigricans",
               "Hypoplectrus_guttavarius", "Hypoplectrus_indigo", "Hypoplectrus_puella",
               "Serranus_tortugarum", "Serranus_tabacarius", "Schultzea_beta",
               "Diplectrum_formosum", "Diplectrum_bivittatum", "Diplectrum_pacificum",
               "Diplectrum_maximum", "Serranus_notospilus", "Serranus_phoebe",
               "Serranus_psittacinus", "Serranus_baldwini", "Serranus_tigrinus",
               "Paralabrax_albomaculatus", "Paralabrax_dewegeri", "Paralabrax_callaensis",
               "Paralabrax_loro", "Paralabrax_auroguttatus", "Paralabrax_clathratus",
               "Paralabrax_humeralis", "Paralabrax_nebulifer", "Paralabrax_maculatofasciatus",
               "Zalanthias_kelloggi", "Serranus_cabrilla", "Serranus_atricauda",
               "Serranus_scriba", "Serranus_hepatus", "Serranus_accraensis",
               "Centropristis_striata", "Chelidoperca_occipitalis", "Chelidoperca_investigatoris",
               "Chelidoperca_pleurospilus")

edvr_serr <- edvr %>% 
  subtreeBAMM(tips = serranids)

edvr_serr_short <- edvr_serr
edvr_serr_short$tip.label <- edvr_serr$tip.label %>% 
  str_replace(pattern = "([A-Z])[a-z]*_([a-z]*)", "italic(\\1.~\\2)")%>% 
  str_replace(pattern = "beta", "'beta'") %>% 
  ggplot2:::parse_safe()

clr_tree <- scico::scico(6, palette = "berlin") %>% # viridis::plasma(6) %>% #
  prismatic::clr_desaturate(shift = .4) %>% 
  # prismatic::clr_lighten(shift = .2) #%>% 
  prismatic::clr_darken(shift = .2)
# plot(edvr_serr_short, labels = TRUE, legend = TRUE,pal = rainbow(15),logcolor = TRUE, cex = .7)

p1 <- as.grob(function(){
  par(mar = c(0,0,0,0))
  bammplot_k(x = edvr_serr_short,
             labels = T,
             lwd = .8, #legend = TRUE,
             cex = .3,
             pal = clr_tree,
             labelcolor = c(rep("black", 9),
                            rep("darkgray", 31)))
  leg_shift_x <- 0
  lines(x = c(0,25) + leg_shift_x,
        y = c(1.5, 1.5),
        col = "darkgray")
  text(x = 12.5 + leg_shift_x,
       y = .5,
       labels = "25 MYR",
       cex = .4, 
       col = "darkgray")
}
)

c1 <- "transparent"
c2 <- rgb(0, 0, 0, .1)
c3 <- rgb(0, 0, 0, .2)
grad_mat <- c(c1, c2, c3, c3) 
dim(grad_mat) <- c(1, length(grad_mat))
grob_grad <- rasterGrob(grad_mat,
                        width = unit(1, "npc"),
                        height = unit(1, "npc"), 
                        interpolate = TRUE)

blank_hamlet <- hypoimg::hypo_outline %>% 
  ggplot()+
  coord_equal()+
  # geom_polygon(aes(x, y), color = "black", fill = rgb(0, 0, 0, .2), size = .1)+
  geom_polygon(aes(x, y), color = rgb(0, 0, 0, .5), fill = rgb(1, 1, 1, .3), size = .1)+
  theme_void()

p_tree <- ggplot() +
  geom_point(data = tibble(v = c(.056, 2.4)),
             x = .5, y = .5, aes(color = v),alpha = 0)+
  scale_color_gradientn(colours = clr_tree, limits = c(.056, 2.4))+
  annotation_custom(grob = grob_grad,
                    ymin = 0, ymax = .22,
                    xmin = .4, xmax = .96)+
  annotation_custom(grob = ggplotGrob(blank_hamlet),
                    xmin = 0.55, xmax = .75,
                    ymin = 0.015, ymax = .14)+
  annotation_custom(grob = p1,
                    xmin = -.16,
                    xmax = 1.05,
                    ymin = -.22,
                    ymax = 1) +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(0, 1),
                  expand = 0)+
  guides(color = guide_colorbar(title = "Diversification Rate",
                                title.position = "top",
                                direction = "horizontal",
                                barheight = unit(3, "pt"),
                                barwidth = unit(57, "pt"),
                                ticks.colour = "white")) +
  # theme_minimal()+
  theme_void()+
  theme(legend.position = c(.01, .08),
        legend.justification = c(0, 0))

# =================================================================================================

plot_fish2 <- function (name, x = 0, y = 3, height = 4, width = 4, lwd = .15, line_color = "black") { 
  spec <- str_remove(string = name, "Hypoplectrus_")
  hypo_anno_l_lwd(spec, xmin = x - 0.5 * width, xmax = x + 
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 * 
                    height, lwd = lwd, line_color = line_color)
}

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

plot_serranus <- function(grob, x = 0, y = 3, height = 4, width = 4, ...){
  annotation_custom(grob = grob, xmin = x - 0.5 * width, xmax = x + 
                      0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 * 
                      height)
}

# import Fst
fst_files <- dir(fst_dir, pattern = '.50k.windowed.weir.fst.gz')

fst_data <- str_c(fst_dir, fst_files) %>%
  purrr::map(summarize_fst) %>%
  bind_rows()

# determine fst ranking
fst_order <- fst_data %>%
  select(run, `mean_weighted-fst`) %>%
  mutate(run = fct_reorder(run, `mean_weighted-fst`))

fst_data_gatger <- fst_data %>% 
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
 
# sort run by average genome wide Fst
run_ord <- tibble(run = levels(fst_data_gatger$run),
                  run_ord = 1:length(levels(fst_data_gatger$run)))

clr_alt <- clr
clr_alt[["uni"]] <- "lightgray"
pca_fish_scale <- 1.15

pca_fish_pos <- tibble(pop = GenomicOriginsScripts::pop_levels,
                       short = str_sub(pop, 1, 3),
                       loc = str_sub(pop, 4, 6),
                       width = c(bel = .08, hon =.08, pan = .09)[loc] * pca_fish_scale,
                       height = c(bel = .08, hon =.08, pan = .09)[loc] * pca_fish_scale) %>% 
  arrange(loc) %>%
  mutate(x = c(-.18, -.01, .03, -.03, .075,
               -.04, -.2, -.02, -.01, -.02, -.01,
               -.15, 0, .05),
         y = c(.02, .27, -.13, -.03, .05,
               -.1, .05, 0, .075, -.23, .2,
               .06, -.2, .2)) %>%
  select(-pop) %>%
  group_by(loc) %>%
  nest()

plot_fish_lwd <- function (short, x = 0, y = 3, height = 5, width = 5, lwd = .15, line_color = "transparent") { 
  hypo_anno_l_lwd(sp_names[short], xmin = x - 0.5 * width, xmax = x + 
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 * 
                    height, lwd = lwd, line_color = line_color)
}

pca_plot <- function(loc){
  set.seed(42)
  evs <- str_c("2_analysis/pca/", loc ,".exp_var.txt.gz") %>% 
    read_tsv()
  str_c("2_analysis/pca/", loc ,".scores.txt.gz") %>% 
    read_tsv() %>% 
    mutate(spec = str_sub(id, -6,-4)) %>%
    ggplot(aes(x = EV01, y = EV02, fill = spec))+
    ggforce::geom_mark_ellipse(aes(color = spec),
    fill = "transparent",
    linetype = 3,
    size = .3,
    expand = unit(5, "pt"))+
    geom_point(shape = 21, aes(color = after_scale(prismatic::clr_darken(fill))), size = .7) +
    (pca_fish_pos$data[[which(pca_fish_pos$loc == loc)]] %>%
       pmap(plot_fish_lwd))+
    labs(x = str_c("PC1 (", sprintf("%.2f",evs$exp_var[[1]]), " %)"),
         y = str_c("PC2 (", sprintf("%.2f",evs$exp_var[[2]]), " %)"))+
    scale_fill_manual(values = clr)+
    scale_color_manual(values = clr_alt %>%
                         prismatic::clr_alpha(alpha = .7) %>%
                         set_names(nm = names(clr_alt)))+
    labs(title = loc_names[[loc]])+
    theme_minimal()+
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          panel.background = element_rect(color = clr_loc[[loc]],
                                          fill = "transparent"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = plot_text_size),
          plot.title = element_text(color = clr_loc[[loc]]))
}

pcas <- c("bel", "hon", "pan") %>% map(pca_plot)
fst_sig_attach <- read_tsv(fst_permutation_file)

# assemble panel b
p2 <- fst_data_gatger %>%
  filter(popstat == "weighted-fst") %>%
  left_join(fst_sig_attach) %>% 
  mutate(loc = str_sub(run, -3, -1)) %>% 
  group_by(loc) %>% 
  mutate(loc_n = length(loc),
         is_sig = p_perm > .05/loc_n) %>% 
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
  geom_point(aes(x = x, y = mean, shape = is_sig, fill = after_scale(color)),
             size = .8)+
  scale_x_continuous(name = "Index of Sympatric Hamlet Species Pair",
                     breaks = 1:28) +
  scale_y_continuous(#breaks = c(0,.05,.1,.15),
    name = expression(italic(F[ST])))+
  scale_color_manual(values = c(make_faint_clr('bel'),
                                make_faint_clr('hon'),
                                make_faint_clr('pan'))[c(2, 4, 6)])+
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 21)) +
  coord_cartesian(xlim = c(0,29),
                  expand = c(0,0))+
  theme_minimal()+
  theme(text = element_text(size = plot_text_size),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y.right = element_text(color = clr_sec),
        axis.title.y.right = element_text(color = clr_sec))

p_done <- (wrap_elements(plot = p_tree +
             theme(axis.title = element_blank(),
                   text = element_text(size = plot_text_size)
                   # plot.margin = unit(c(0, 0, 0,0), "pt"),
                   # plot.background = element_rect(fill = "green")
                   ),
             clip = FALSE) +
             p2) /
  (pcas %>% wrap_plots()) +
  plot_layout(heights = c(1,.75)) +
  plot_annotation(tag_levels = 'a') &
  theme(text = element_text(size = plot_text_size),
        plot.background = element_rect(fill = "transparent",
                                       color = "transparent"))

scl <- .75
hypo_save(p_done, filename = 'figures/F1_redesign.pdf',
          width = 9 * scl,
          height = 6 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
