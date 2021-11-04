#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFxx3.R 
# ===============================================================
# This script produces Figure 4 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/astral/astral_5000x_5kb_v1_all.tre")
# script_name <- "R/fig/plot_SFxx3.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)
library(tidygraph)
library(ggraph)
library(patchwork)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
# config -----------------------
tree_hypo_file <- as.character(args[1])
tree <- read.tree(tree_hypo_file)
tree$edge.length <- replace(tree$edge.length, tree$edge.length == "NaN", 0.05)   # Set terminal branches to 0.05
tree$edge.length[c( 81, 83)] <- tree$edge.length[c( 81, 83)] * 0.1
# which(tree$edge.length > 1)



tree_rooted <- root(phy = tree, outgroup = c("s_tort_3torpan", "20478tabhon", "28393torpan"))
# tree_s_mid <- phangorn::midpoint(tree_rooted)
clr_neutral <- rgb(.2, .2, .2)

### Prepare tree and categorize support values
tree_plus <- ggtree(tree_rooted)  %>%
  .$data %>%
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
         support = as.numeric(label) * 100,
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
          `branch.length` = if_else(node %in% c( 212, 213), `branch.length` * .0001, `branch.length`),
          branch_type = if_else(node %in% c(212, 213), "broken", "whole")
         )
  
tree_plus %>%  arrange(branch_type, -branch.length)
ggtree(tree_plus)

(t_plot <- (ggtree(tr = tree_plus,
                   layout = 'fan', #open.angle = 180,
                      aes(color = spec, linetype = branch_type), size = .2)) + #%>% 
       # ggtree::rotate_tree(angle = -100)) +
    geom_tippoint(aes(color = spec,
                      shape = loc,
                      fill = after_scale(color)), size = .5) +
    geom_nodepoint(data = tree_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = clr_neutral) +
    scale_color_manual(values = c(GenomicOriginsScripts::clr2, ungrouped = "gray60"), labels = GenomicOriginsScripts::sp_labs
                       #guide = 'none'
                       ) +
    scale_shape_manual(values = c(bel = 21, flo = 24, hon = 22, pan = 23), labels = GenomicOriginsScripts::loc_names#,
                       #guide = 'none'
                       ) +
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
    scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
    # Add scale bar:
    ggtree::geom_treescale(width = .2,
                           x = .13, y = 85.5,
                           offset = -7,
                           linesize = .2,
                           fontsize = plot_text_size/.pt,
                           color = clr_neutral) +
    # scale_x_continuous(limits = c(-.05, .26), expand = c(0,0)) +
    # xlim(c(-.05,#-.15,
    #        .265)) +
    guides(fill = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2, label.hjust = 0),
           size = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2, label.hjust = 0),
           shape = guide_legend(title = "Location", title.position = "top",
                               nrow = 2, label.hjust = 0),
           color = guide_legend(title = "Species", title.position = "top",
                               ncol = 2, label.hjust = 0)) +
    theme_void() +
    theme(legend.position = 'bottom',
      legend.title.align = 0,
      legend.text = element_text(color = "gray20"),
      legend.title = element_text(color = "gray20")) )

# clado_plus <- ggtree(tree_s_mid, branch.length = 'none',
#                      layout = 'fan', open.angle = 180)  %>%
#   .$data %>%
#   mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
#          loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
#          support = as.numeric(label) * 100,
#          support_class = cut(support, c(0,50,70,90,100)) %>% 
#            as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
#   ) 
# 
# (c_plot <- ggtree(tr = clado_plus,
#        layout = 'fan', #open.angle = 180,
#        aes(color = spec), size = .2)+
#   geom_tippoint(aes(color = spec,
#                     shape = loc,
#                     fill = after_scale(color)), size = .5) +
#   geom_nodepoint(data = clado_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
#                  aes(fill = support_class,
#                      size = support_class),
#                  shape = 21,
#                  color = clr_neutral) +
#   scale_color_manual(values = c(GenomicOriginsScripts::clr2, ungrouped = "gray60"),
#                      guide = 'none') +
#   scale_shape_manual(values = c(bel = 21, flo = 24, hon = 22, pan = 23), labels = GenomicOriginsScripts::loc_names,
#                      guide = 'none') +
#   scale_fill_manual(values = c(`(0,50]`   = "transparent",
#                                `(50,70]`  = "white",
#                                `(70,90]`  = "gray",
#                                `(90,100]` = "black"),
#                     drop = FALSE) +
#   scale_size_manual(values = c(`(0,50]`   = 0,
#                                `(50,70]`  = .8,
#                                `(70,90]`  = .8,
#                                `(90,100]` = .8),
#                     na.value = 0,
#                     drop = FALSE) +
#   # Add scale bar:
#   # ggtree::geom_treescale(width = 1,
#   #                        x = .13, y = 170,
#   #                        offset = -7,
#   #                        linesize = .2,
#   #                        fontsize = plot_text_size/.pt,
#   #                        color = clr_neutral) +
#   # scale_x_continuous(limits = c(-.05, .26), expand = c(0,0)) +
#   # xlim(c(-.05,#-.15,
#   #        .265)) +
#   guides(fill = guide_legend(title = "Node Support Class", title.position = "top",
#                              nrow = 1),
#          size = guide_legend(title = "Node Support Class", title.position = "top",
#                              nrow = 1)) +
#   theme_void() +
#   theme(#legend.position = 'bottom',
#     legend.title.align = 0.5,
#     legend.text = element_text(color = "gray20"),
#     legend.title = element_text(color = "gray20")) )

y_sep <- .5
x_shift <- .85

( p_done <- ggplot() +
    coord_equal(xlim = c(0, 1),
                ylim = c(0, 1),
                expand = 0) +
    annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
                      ymin = -1 - y_sep , 
                      ymax = .7 + y_sep,
                      xmin = 0 - x_shift,
                      xmax = 1 + x_shift) +
    annotation_custom(grob = cowplot::get_legend(t_plot),
                      ymin = -.1, ymax = .3,
                      xmin = 0, xmax = 1) +
    
    # theme_dark()
    theme_void()
)

# p_done <- ggplot() +
#     coord_equal(xlim = c(0, 1),
#                 ylim = c(0, 1),
#                 expand = 0) +
#     annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
#                       ymin = -.35 - y_sep , 
#                       ymax = .7 + y_sep,
#                       xmin = 0 - x_shift,
#                       xmax = 1 + x_shift) +
#     annotation_custom(grob = cowplot::get_legend(t_plot),
#                       ymin = -.1, ymax = .3,
#                       xmin = 0, xmax = 1) +
#     theme_void()

y_sep <- .1
x_shift <- .1
(p_tdone <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
                    ymin = 0 - y_sep , 
                    ymax = 1 + y_sep,
                    xmin = 0 - x_shift,
                    xmax = 1 + x_shift) +
  theme_void())

p_cdone <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(c_plot + theme(legend.position = "none")),
                    ymin = 0 - y_sep , 
                    ymax = 1 + y_sep,
                    xmin = 0 - x_shift,
                    xmax = 1 + x_shift) +
  theme_void()

# p_done <- cowplot::plot_grid(cowplot::plot_grid(p_tdone, p_cdone,labels = letters[1:2],label_size = plot_text_size, label_fontface = "plain"),
#                              cowplot::get_legend(t_plot +
#                                                    theme_minimal(base_size = plot_text_size) + 
#                                                    theme(legend.position = "bottom",
#                                                          legend.title.align =0,
#                                                          legend.key.height = unit(7,"pt"))),
#                              ncol = 1, rel_heights = c(1,.2))

p_done <- cowplot::plot_grid(p_tdone, cowplot::get_legend(t_plot +
                                                  theme_minimal(base_size = plot_text_size) + 
                                                  theme(legend.position = "right",
                                                        legend.title.align =0,
                                                        legend.key.height = unit(7,"pt"))),rel_widths = c(1,.5))
scl <- .75
hypo_save(p_done, filename = 'figures/SFxx3_2.pdf',
          width = f_width,
          height = .6 * f_width,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
