#!/usr/bin/env Rscript
#
# Context: this script depends on the input file 2_analysis/summaries/fst_permutation_summary.tsv
#          which is created by R/fst_permutation.R
#
# run from terminal:
# Rscript --vanilla R/fig/plot_F1.R \
#   2_analysis/fst/50k/ \
#   2_analysis/summaries/fst_globals.txt \
#   2_analysis/summaries/fst_permutation_summary.tsv \
#   2_analysis/fotl/concat_R24ed.treefile
# ===============================================================
# This script produces Figure 1 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt',
#           '2_analysis/summaries/fst_permutation_summary.tsv',
#           "2_analysis/fotl/concat_R24ed.treefile")
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
library(ggtree)

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
tree_file <- as.character(args[4])
wdh <- .3          # The width of the boxplots
scaler <- 20       # the ratio of the Fst and the dxy axis (legacy - not really needed anymore)
clr_sec <- 'gray'  # the color of the secondary axis (dxy)

# start script -------------------
tree <- read.tree(tree_file) 
tree_rooted <- root(phy = tree, outgroup = "Epinephelus_maculatus")
clr_neutral <- rgb(.2, .2, .2)
clr_highlight <- "gray40" #"#FF8029"
clr_stars <- "firebrick" 

### Edit tip labels
tree_rooted$tip.label <- tree_rooted$tip.label %>% 
  str_replace(pattern = "20864abehon", "Hypoplectrus_aberrans") %>%
  str_replace(pattern = "20642gumhon", "Hypoplectrus_gummigutta") %>%
  str_replace(pattern = "18238indbel", "Hypoplectrus_indigo") %>%
  str_replace(pattern = "PL17_122maybel", "Hypoplectrus_maya") %>%
  str_replace(pattern = "18906nigpan", "Hypoplectrus_nigricans") %>%
  str_replace(pattern = "18434puepan", "Hypoplectrus_puella") %>% 
  str_replace(pattern = "20613ranhon", "Hypoplectrus_randallorum") %>%
  str_replace(pattern = "18448unipan", "Hypoplectrus_unicolor") %>%
  str_replace(pattern = "PL17_160floflo", "Hypoplectrus_floridae") %>%
  str_replace(pattern = "20478tabhon", "Serranus_tabacarius") %>%
  str_replace(pattern = "s_tort_3torpan", "Serranus_tortugarum") %>% 
  #
  str_replace(pattern = "([A-Z])([a-z])[a-z]*_([a-z]*)", "\\1\\2. \\3")  %>%
  str_replace(pattern = "Ce.", "Cp.")%>%
  str_replace(pattern = "Za.", "Pl.") %>%
  str_replace(pattern = "Hy.", "H.") %>% 
  str_replace(pattern = "Di.", "D.") %>% 
  str_replace(pattern = "Ep.", "E.")

### Prepare tree, categorize support values and define group
tree_plus <- ggtree(tree_rooted, layout = "rectangular", ladderize = TRUE, right = TRUE) %>%
  #flip(25,26) %>%
  .$data %>%
  mutate(support = as.numeric(label),
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")))

hamlets <- tree_plus$label[grepl(pattern = "H. ", tree_plus$label)]

tree_plus <- tree_plus %>% 
  groupOTU(.node = hamlets)

### Define groups
our_taxa <- c("H. aberrans", "H. gummigutta", "H. indigo", "H. maya", "H. nigricans", "H. puella", 
              "H. randallorum", "H. unicolor", "H. floridae", "Se. tabacarius", "Se. tortugarum")

hermaphrodites <- tree_plus %>% filter(node %in% c(60, 63)) %>% 
  mutate(x = if_else(node == 60, x - branch.length, x - branch.length *.5),
         y = if_else(node == 60, .5 * (y + tree_plus$y[tree_plus$node == 61]), y),
         star = "\U2605")

# c1 <- "transparent"
c2 <- prismatic::clr_alpha(clr_highlight, .1)
# c3 <- prismatic::clr_alpha(clr_highlight, .2)

# grad_mat <- c(c1, c2, c3, c3) 
# dim(grad_mat) <- c(1, length(grad_mat))
# grob_grad <- rasterGrob(grad_mat,
#                         width = unit(1, "npc"),
#                         height = unit(1, "npc"), 
#                         interpolate = TRUE)

grob_rect <- rectGrob(gp = grid::gpar(fill = c2, col = "transparent"))

blank_hamlet <- hypoimg::hypo_outline %>% 
  ggplot()+
  coord_equal()+
  geom_polygon(aes(x, y), color = rgb(0,0,0,.3), fill = rgb(1, 1, 1, .3), size = .1)+
  theme_void()

### Draw tree
(p_tree <- ggtree(tree_plus,
                 aes(color = group), size = .2) +
  annotation_custom(grob = grob_rect,#grad,
                    ymin = .2, ymax = 12.5,
                    xmin = .042, xmax = .125)+
  annotation_custom(grob = ggplotGrob(blank_hamlet),
                    xmin = 0.09, xmax = .125,
                    ymin = 1.2, ymax = 10.5) +
  geom_tiplab(aes(color = group, label = if_else(label %in% our_taxa, str_c(label,"*"), label)),   # add asterisks to our taxa
              size = 1.3, hjust = -.1,
              fontface = 'italic') +
  ggplot2::xlim(0, 0.12) +   # add extra space for long labels
  geom_nodepoint(data = tree_plus %>% filter(!is.na(support_class)), 
                   aes(fill = support_class,
                     size = support_class),
                 shape = 21, color = clr_neutral) +
  geom_text(data = hermaphrodites, aes(x = x, y = y, label = star),
            family = "DejaVu Sans", color = clr_stars,
            size = 2, vjust = .35) +
  scale_color_manual(values = c(clr_neutral, clr_neutral, clr_highlight)) +
  scale_fill_manual(values = c(`(0,50]`   = "transparent",
                               `(50,70]`  = "white",
                               `(70,90]`  = "gray",
                               `(90,100]` = "black"),
                    drop = FALSE) +
  scale_size_manual(values = c(`(0,50]`   = 0,
                               `(50,70]`  = .8,
                               `(70,90]`  = .8,
                               `(90,100]` = .8),
                    na.value = 0,
                    drop = FALSE) +
  geom_treescale(color = clr_neutral, 
                 fontsize = 2, linesize = .2, 
                 x = 0.045, y = 1.5) +
  guides(fill = guide_legend(title = "Node Support Class", title.position = "top", override.aes = list(color = clr_neutral), nrow = 2),
         size = guide_legend(title = "Node Support Class", title.position = "top", override.aes = list(color = clr_neutral), nrow = 2),
         color = 'none') +
    coord_cartesian(xlim = c(-.005,.125), ylim = c(0,45), expand = 0) +
    theme_void() +
  theme(legend.position = c(0.05,0),
        legend.justification = c(0,0),
        legend.title.align = 0,
        legend.key.height = unit(8,"pt"),
        legend.key.width = unit(6,"pt"),
        legend.text = element_text(color = clr_neutral),
        legend.title = element_text(color = clr_neutral))
)

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
