#!/usr/bin/env Rscript
#
# Context: this script depends on the input file 2_analysis/summaries/fst_permutation_summary.tsv
#          which is created by R/fst_permutation.R
#
# run from terminal:
# Rscript --vanilla R/fig/plot_F1.R \
#     2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt 2_analysis/summaries/fst_permutation_summary.tsv
# ===============================================================
# This script produces Figure 1 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt',
#           '2_analysis/summaries/fst_permutation_summary.tsv')
# script_name <- "R/fig/plot_F1.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
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

label_two_chars <- c(`italic(Z.~kelloggi)` = "italic(Pl.~kelloggi)",
                     `italic(P.~maculatofasciatus)` = "italic(Pa.~maculatofasciatus)",
                     `italic(P.~nebulifer)` = "italic(Pa.~nebulifer)",
                     `italic(P.~humeralis)` = "italic(Pa.~humeralis)",
                     `italic(P.~clathratus)` = "italic(Pa.~clathratus)",
                     `italic(P.~auroguttatus)` = "italic(Pa.~auroguttatus)",
                     `italic(P.~loro)` = "italic(Pa.~loro)",
                     `italic(P.~callaensis)` = "italic(Pa.~callaensis)",
                     `italic(P.~dewegeri)` = "italic(Pa.~dewegeri)",
                     `italic(P.~albomaculatus)` = "italic(Pa.~albomaculatus)",
                     `italic(S.~'beta')` = "italic(Sc.~'beta')",
                     `italic(S.~notospilus)` = "italic(Se.~notospilus)",
                     `italic(S.~phoebe)` = "italic(Se.~phoebe)",
                     `italic(S.~psittacinus)` = "italic(Se.~psittacinus)",
                     `italic(S.~baldwini)` = "italic(Se.~baldwini)",
                     `italic(S.~tigrinus)` = "italic(Se.~tigrinus)",
                     `italic(S.~cabrilla)` = "italic(Se.~cabrilla)",
                     `italic(S.~atricauda)` = "italic(Se.~atricauda)",
                     `italic(S.~scriba)` = "italic(Se.~scriba)",
                     `italic(S.~hepatus)` = "italic(Se.~hepatus)",
                     `italic(S.~accraensis)` = "italic(Se.~accraensis)",
                     `italic(C.~striata)` = "italic(Cp.~striata)",
                     `italic(C.~occipitalis)` = "italic(Ch.~occipitalis)",
                     `italic(C.~investigatoris)` = "italic(Ch.~investigatoris)",
                     `italic(C.~pleurospilus)` = "italic(Ch.~pleurospilus)",
                     `italic(S.~tabacarius)` = "italic(Se.~tabacarius)",
                     `italic(S.~tortugarum)` = "italic(Se.~tortugarum)")

edvr_serr_short <- edvr_serr
edvr_serr_short$tip.label <- edvr_serr$tip.label %>% 
  str_replace(pattern = "([A-Z])[a-z]*_([a-z]*)", "italic(\\1.~\\2)")%>% 
  str_replace(pattern = "beta", "'beta'") %>%
  ifelse(. %in% names(label_two_chars), label_two_chars[.], .) %>% 
  ggplot2:::parse_safe()

clr_tree <- scico::scico(6, palette = "berlin") %>% 
  prismatic::clr_desaturate(shift = .4) %>% 
  prismatic::clr_darken(shift = .2)


clr_lab <- rep(c("black","darkgray"), c(9,31))
# clr_lab <- rep(c("black","darkgray","black","darkgray"), c(21,10,4,5))

p1 <- as.grob(function(){
  par(mar = c(0,0,0,0))
  bammplot_k(x = edvr_serr_short,
             labels = T,
             lwd = .8,
             cex = .3,
             pal = clr_tree,
             labelcolor = clr_lab)
  leg_shift_x <- 0
  text(x = c(21.2, 40.2), y = c(15.6, 33.25),
       label = "\U2605", family = "DejaVu Sans", col = clr_tree[[6]], cex = .5)
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
  geom_polygon(aes(x, y), color = rgb(0, 0, 0, .5), fill = rgb(1, 1, 1, .3), size = .1)+
  theme_void()

p_tree <- ggplot() +
  geom_point(data = tibble(v = c(.056, 2.4)),
             x = .5, y = .5, aes(color = v),alpha = 0)+
  scale_color_gradientn(colours = clr_tree, limits = c(.056, 2.4))+
  annotation_custom(grob = grob_grad,
                    ymin = 0.01, ymax = .2335,
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
  guides(color = guide_colorbar(title = "Speciation Rate",
                                title.position = "top",
                                direction = "horizontal",
                                barheight = unit(3, "pt"),
                                barwidth = unit(57, "pt"),
                                ticks.colour = "white")) +
  theme_void()+
  theme(legend.position = c(.01, .08),
        legend.justification = c(0, 0))

globals <- vroom::vroom(fst_globals, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  mutate(run = str_c(str_sub(run,1,3),loc,'-',str_sub(run,5,7),loc),
         run = fct_reorder(run,weighted))

# sort run by average genome wide Fst
run_ord <- tibble(run = levels(globals$run),
                  run_ord = seq_along(levels(globals$run)))

# load fst permutation results
fst_sig_attach <- read_tsv(fst_permutation_file) %>% 
  filter( subset_type == "whg" ) %>% 
  mutate(loc = str_sub(run, -3, -1)) %>%
  group_by(loc) %>% 
  mutate(loc_n = 28,#length(loc),
         fdr_correction_factor =  sum(1 / 1:loc_n),
         fdr_alpha = .05 / fdr_correction_factor,
         is_sig = p_perm > fdr_alpha) %>% 
  ungroup()

# create network annotation
# underlying structure for the network plots
networx <- tibble( loc = c('bel','hon', 'pan'),
                   n = c(5, 6, 3),
                   label = list(str_c(c('ind','may','nig','pue','uni'),'bel'),
                                str_c(c('abe','gum','nig','pue','ran','uni'),'hon'),
                                str_c(c('nig','pue','uni'),'pan')),
                   weight = c(1,1.45,1)) %>%
  purrr::pmap_dfr(network_layout) %>%
  mutate(edges = map(edges, function(x){x %>% left_join(globals,by = "run") }))

# plot the individual networks by location
plot_list <- networx %>%
  purrr::pmap(plot_network, node_lab_shift = .2)

# combine the networks into a single grob
p_net <- cowplot::plot_grid(
  plot_list[[1]] + theme(legend.position = "none"),
  plot_list[[2]] + theme(legend.position = "none"),
  plot_list[[3]] + theme(legend.position = "none"),
  ncol = 3) %>%
  cowplot::as_grob()

p2 <- globals %>%
  left_join(fst_sig_attach %>% mutate(run = factor(run, levels = levels(globals$run))) ) %>% 
  ggplot(aes(color = loc)) +
  geom_bar(aes(x = as.numeric(run),
                   y = weighted,
               alpha = is_sig,
               fill = after_scale(clr_lighten(color))),
               stat = "identity",size = .2, width = .8)+
  annotation_custom(p_net, ymin = .05, xmax = 24.5) +
  scale_x_continuous(name = "Pair of sympatric species",
                     breaks = 1:28) +
  scale_y_continuous(name = expression(italic(F[ST])))+
  scale_color_manual(values = c(make_faint_clr('bel'),
                                make_faint_clr('hon'),
                                make_faint_clr('pan'))[c(2, 4, 6)])+
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 21)) +
  scale_alpha_manual(values = c(`TRUE` = .15, `FALSE` = 1)) + 
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
# assemble panel c-e
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

pcas <- c("bel", "hon", "pan") %>% map(pca_plot)

fish_tib <- tibble(short = names(clr)[!names(clr) %in% c("flo", "tab", "tor")],
       x = c(0.5,  3.5,  7,  9.7, 12.25, 15.25, 18, 21.5)
       )
# ===================
key_sz <- .75
p_leg <- fish_tib %>% 
  ggplot() +
  coord_equal(xlim = c(-.05, 24), expand = 0) +
  geom_tile(aes(x = x, y = 0,
                fill = short, 
                color = after_scale(prismatic::clr_darken(fill, .25))),
            width = key_sz, height = key_sz, size = .3) +
  geom_text(aes(x = x + .6, y = 0,
                label = str_c("H. ", sp_names[short])), 
            hjust = 0, fontface = "italic", size = plot_text_size / ggplot2:::.pt) +
  pmap(fish_tib, plot_fish_lwd, width = 1, height = 1, y = 0) +
  scale_fill_manual(values = clr, guide = FALSE) +
  theme_void()


p_combined <- ((wrap_elements(plot = p_tree +
                                theme(axis.title = element_blank(),
                                      text = element_text(size = plot_text_size)),
                              clip = FALSE) + p2 ) /
                 (pcas %>% wrap_plots()) +
                 plot_layout(heights = c(1,.75)) +
                 plot_annotation(tag_levels = 'a') &
                 theme(text = element_text(size = plot_text_size),
                       plot.background = element_rect(fill = "transparent",
                                                      color = "transparent"))) 

p_done <- cowplot::plot_grid(p_combined, p_leg, ncol = 1, rel_heights = c(1,.06))

scl <- .75
hypo_save(p_done, filename = 'figures/F1.pdf',
          width = 9 * scl,
          height = 6 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
