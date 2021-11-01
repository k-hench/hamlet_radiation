#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F4.R \
#     2_analysis/raxml/hyp155_n_0.33_mac4_5kb.raxml.support \
#     2_analysis/ibd/no_outgr_direct_10.ibd.tsv
# ===============================================================
# This script produces Figure 4 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
args <- c("2_analysis/astral/astral_5000x_5kb_v1_noS.tre",
          "2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv")
# script_name <- "R/fig/plot_F4.R"
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
ibd_file <- as.character(args[2])

tree <- read.tree(tree_hypo_file)
tree$edge.length <- replace(tree$edge.length, tree$edge.length == "NaN", 0.05)   # Set terminal branches to 0.05

tree_rooted <- root(phy = tree, outgroup = "PL17_160floflo")
clr_neutral <- rgb(.2, .2, .2)


### Prepare tree and categorize support values
tree_plus <- ggtree(tree_rooted) %>%
  # flip(77,78) %>%
  # flip(196,205) %>%
  # ggtree::rotate(node = 195) %>% 
  ggtree::rotate(node = 205) %>% 
  # ggtree::rotate(node = 227) %>% 
  # ggtree::rotate(node = 228) %>% 
  .$data %>%
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
         support = as.numeric(label) * 100,
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
  ) 


# (t_plot <- (ggtree(tree_plus, layout = 'fan', open.angle = 180,
#                   aes(color = spec), size = .2)) + #%>% ggtree::rotate_tree(angle = -30))+
#   # geom_tiplab2(aes(color = lab2spec(label),
#   #                  label = str_sub(
#   #                  label, -6, -1)),
#   #              size = 3, hjust = -.1) +
#   geom_tippoint(aes(color = spec,
#                     shape = loc,
#                     fill = after_scale(color)), size = .5) +
#   geom_nodepoint(data = tree_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
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
#   # ggtree::geom_treescale(width = .002,
#   #                        x = -.0007, y = 155, 
#   #                        offset = -3,fontsize = 3,
#   #                        color = clr_neutral) +
#   xlim(c(0,#-.05,#-.15,
#          .265)) +
#   guides(fill = guide_legend(title = "Node Support Class", title.position = "top", nrow = 2),
#          size = guide_legend(title = "Node Support Class", title.position = "top", nrow = 2)) +
#   theme_void() +
#   theme(#legend.position = 'bottom',
#         legend.title.align = 0.5,
#         legend.text = element_text(color = "gray20"),
#         legend.title = element_text(color = "gray20")) )

(t_plot <-
    (ggtree(tree_plus, layout = 'fan', #open.angle = 180,
                   aes(color = spec), size = .2) %>% 
    ggtree::rotate_tree(angle = -100)) +
    geom_tippoint(aes(color = spec,
                      shape = loc,
                      fill = after_scale(color)), size = .5) +
    geom_nodepoint(data = tree_plus %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
                   aes(fill = support_class,
                       size = support_class),
                   shape = 21,
                   color = clr_neutral) +
    scale_color_manual(values = c(GenomicOriginsScripts::clr2, ungrouped = "gray60"),
                       guide = 'none') +
    scale_shape_manual(values = c(bel = 21, flo = 24, hon = 22, pan = 23), labels = GenomicOriginsScripts::loc_names,
                       guide = 'none') +
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
    # Add scale bar:
    ggtree::geom_treescale(width = .05,
                           x = .1, y = 130,
                           offset = -1,
                           linesize = .2,
                           fontsize = plot_text_size/.pt,
                           color = clr_neutral) +
    scale_x_continuous(limits = c(-.05, .26), expand = c(0,0)) +
    # xlim(c(-.05,#-.15,
    #        .265)) +
    guides(fill = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2),
           size = guide_legend(title = "Node Support Class", title.position = "top",
                               nrow = 2)) +
    theme_void() +
    theme(#legend.position = 'bottom',
      legend.title.align = 0.5,
      legend.text = element_text(color = "gray20"),
      legend.title = element_text(color = "gray20")) 
  )

y_sep <- .55
x_shift <- .5

( p1 <- ggplot() +
  coord_equal(xlim = c(0, 1),
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
                    ymin = -.15 - y_sep , 
                    ymax = 1 + y_sep,
                    xmin = .25 - x_shift,
                    xmax = 1 + x_shift) +
  # annotation_custom(grob = cowplot::get_legend(t_plot),
  #                   ymin = 0.12, ymax = .22,
  #                   xmin = 0, xmax = .53) +

    # theme_dark()
  theme_void()
  )

# p1 <- ggplot() +
#     coord_equal(xlim = c(0.45, 0.46),
#                 ylim = c(0, 0.47),
#                 expand = 0) +
#     annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
#                       ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
#                       xmin = -.1,
#                       xmax = 1.1) +
#     # annotation_custom(grob = cowplot::get_legend(t_plot),
#     #                   ymin = 0.12, ymax = .22,
#     #                   xmin = 0, xmax = .53) +
#     # theme_dark()
#   theme_void()

# (t_final <- ggplot() +
#     coord_equal(xlim = c(0.38, 0.46),
#                 ylim = c(-0.1, 0.42),
#                 expand = 0) +
#     annotation_custom(grob = ggplotGrob(t_plot + theme(legend.position = "none")),
#                       ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
#                       xmin = -.1,
#                       xmax = 1.1) +
#     annotation_custom(grob = cowplot::get_legend(t_plot),
#                       ymin = -0.12, ymax = 0,
#                       xmin = 0, xmax = .93) +
#     theme_dark()
#     # theme_void()
# )

# raxml_tree <- read.tree(tree_hypo_file) 
# raxml_tree_rooted <- root(phy = raxml_tree, outgroup = "PL17_160floflo")
# clr_neutral <- rgb(.6, .6, .6)
# lyout <- 'circular'
# 
# raxml_tree_rooted_grouped <- groupClade(raxml_tree_rooted,
#                                         .node = c(298, 302, 187, 179, 171, 159,
#                                                   193, 204, 201, 222, 219, 209,
#                                                   284, 278, 268, 230, 242),
#                                         group_name =  "clade")
# 
# clade2spec <- c( `0` = "none", `1` = "ran", `2` = "uni", `3` = "ran", `4` = "may",
#                  `5` = "pue", `6` = "ind", `7` = "nig", `8` = "nig", `9` = "ran",
#                  `10` = "abe", `11` = "abe", `12` = "gum", `13` = "uni", `14` = "pue",
#                  `15` = "uni", `16` = "pue", `17` = "nig")
# 
# raxml_data <- ggtree(raxml_tree_rooted_grouped, layout = lyout) %>%
#   .$data %>% 
#   mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
#          support = as.numeric(label),
#          support_class = cut(support, c(0,50,70,90,100)) %>% 
#            as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
#            )
# 
# p_tree <- (open_tree(
#   ggtree(raxml_data, layout = lyout,
#          aes(color = ifelse(clade == 0,
#                             lab2spec(label),
#                             clade2spec[as.character(clade)])), size = .25) %>%
#     ggtree::rotate(200), 180))  +
#   # geom_tippoint(size = .2) + 
#   geom_tiplab2(aes(color = lab2spec(label), label = str_sub(label, -6, -1)),
#   size = GenomicOriginsScripts::plot_text_size_small / ggplot2:::.pt  *.6,#2.5, 
#   hjust = -.1)+
#   ggtree::geom_treescale(width = .002,
#                          linesize = .2,
#                          x = -.0007, y = 155, 
#                          offset = -4,
#                          fontsize = GenomicOriginsScripts::plot_text_size_small / ggplot2:::.pt,
#                          color = clr_neutral) +
#   xlim(c(-.0007,.0092)) +
#   ggtree::geom_nodepoint(aes(fill = support_class, 
#                              size = support_class),
#                  shape = 21#, linewidth = 3
#                  ) +
#   scale_color_manual(values = c(ungrouped = clr_neutral, 
#                                 GenomicOriginsScripts::clr2),
#                      guide = FALSE) +
#   scale_fill_manual(values = c(`(0,50]` = "transparent",
#                                `(50,70]` = "white",
#                                `(70,90]` = "gray",
#                                `(90,100]` = "black"),
#                     drop = FALSE) +
#   scale_size_manual(values = c(`(0,50]` = 0,
#                                `(50,70]` = .4,
#                                `(70,90]` = .4,
#                                `(90,100]` = .4),
#                     na.value = 0,
#                     drop = FALSE)+
#   guides(fill = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2,keyheight = unit(9,"pt")),
#          size = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2,keyheight = unit(9,"pt"))) +
#   theme_void(base_size = GenomicOriginsScripts::plot_text_size_small  ) 
# 
# y_sep <- .05
# x_shift <- -.03
# p1 <- ggplot() +
#   coord_equal(xlim = c(0, .93),
#               ylim = c(-.01, .54),
#               expand = 0) +
#   annotation_custom(grob = ggplotGrob(p_tree + theme(legend.position = "none")),
#                     ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
#                     xmin = -.1, xmax = 1.1) +
#   annotation_custom(grob = cowplot::get_legend(p_tree),
#                     ymin = .35, ymax = .54,
#                     xmin = 0, xmax = .2) +
#   theme_void()

data_ibd <- read_tsv(ibd_file) %>% 
  mutate(ibd_total = (ibd2_cM_m1 + 0.5*ibd1_cM_m1) / (ibd0_cM_m1 + ibd1_cM_m1 + ibd2_cM_m1)) 

set.seed(42)
p2 <- data_ibd %>% 
  as_tbl_graph() %>%
  mutate(spec = str_sub(name, -6, -4),
         loc = factor(str_sub(name, -3, -1), levels = c("bel", "hon", "pan", "flo")))  %>% 
  ggraph( layout = 'fr', weights = ibd_total) +
  geom_edge_link(aes(alpha = ibd_total, edge_width = ibd_total), color = rgb(.1,.1,.1)) +
  geom_node_point(aes(fill = spec,
                      shape = loc, color = after_scale(clr_darken(fill,.3))), size = 1.2) +
  scale_fill_manual("Species", values = GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))],
                    labels = GenomicOriginsScripts::sp_labs)+
  scale_edge_alpha_continuous(#limits = c(0,.1),
    range = c(0,1), guide = "none") +
  scale_edge_width_continuous(#limits = c(0,.1),
    range = c(.1, .4), guide = "none") +
  scale_shape_manual("Site", values = c(bel = 21, flo = 24, hon = 22, pan = 23),
                     labels = GenomicOriginsScripts::loc_names, drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 21, size = 2.5), title.position = "top"),
         shape = guide_legend(nrow = 2, title.position = "top")) +
  coord_equal()  +
  theme(text = element_text(size = GenomicOriginsScripts::plot_text_size),
        panel.background = element_blank())

p_done <- cowplot::plot_grid(cowplot::plot_grid(p1, p2 + theme(legend.position = "none"),rel_widths = c(1,.9),
                                                labels = c("a", "b"), label_fontface = "plain", label_size = plot_text_size),
                             cowplot::plot_grid(cowplot::get_legend(t_plot+ theme_minimal(base_size = GenomicOriginsScripts::plot_text_size)),
                                                cowplot::get_legend(p2 + theme_minimal(base_size = GenomicOriginsScripts::plot_text_size) + theme(legend.position = "bottom")),
                                                rel_widths = c(.3,1)),
                             rel_heights = c(1,.15),
                             ncol = 1)

# p_done <- (p1  + p2 + plot_layout(widths = c(1,1)))  +
#   plot_annotation(tag_levels = "a") +
#   plot_layout(heights = c(1, .02),
#               guides = "collect") &
#   theme(text = element_text(size = GenomicOriginsScripts::plot_text_size),
#         plot.tag.position = c(0, 1),
#         legend.position = "bottom",
#         legend.key = element_blank(),
#         legend.direction = "horizontal",
#         legend.margin = margin(),
#         legend.background = element_blank(),
#         legend.box = "horizontal", 
#         legend.text.align = 0,panel.background = element_rect(fill = "red"))

hypo_save(plot = p_done,
          filename = "figures/Fxx1_2.png",
          width = GenomicOriginsScripts::f_width,
          height = GenomicOriginsScripts::f_width * .65,
          bg = "transparent",
          type = "cairo",
          dpi = 600,
          comment = plot_comment)

system("convert figures/Fxx1.png figures/Fxx1.pdf")
system("rm figures/Fxx1.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/Fxx1.pdf")
system(create_metadata)
