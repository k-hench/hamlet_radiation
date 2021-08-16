#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF9.R 2_analysis/raxml/hyp155_n_0.33_mac4_5kb.raxml.support 2_analysis/ibd/no_outgr_direct_10.ibd.tsv
# ===============================================================
# This script produces Suppl. Figure 9 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/raxml/hyp155_n_0.33_mac4_5kb.raxml.support",
#           "2_analysis/ibd/no_outgr_direct_10.ibd.tsv")
# script_name <- "R/fig/plot_SF9.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
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

raxml_tree <- read.tree(tree_hypo_file) 
raxml_tree_rooted <- root(phy = raxml_tree, outgroup = "PL17_160floflo")
clr_neutral <- rgb(.6, .6, .6)
lyout <- 'circular'

raxml_tree_rooted_grouped <- groupClade(raxml_tree_rooted,
                                        .node = c(298, 302, 187, 179, 171, 159,
                                                  193, 204, 201, 222, 219, 209,
                                                  284, 278, 268, 230, 242),
                                        group_name =  "clade")

clade2spec <- c( `0` = "none", `1` = "ran", `2` = "uni", `3` = "ran", `4` = "may",
                 `5` = "pue", `6` = "ind", `7` = "nig", `8` = "nig", `9` = "ran",
                 `10` = "abe", `11` = "abe", `12` = "gum", `13` = "uni", `14` = "pue",
                 `15` = "uni", `16` = "pue", `17` = "nig")

raxml_data <- ggtree(raxml_tree_rooted_grouped, layout = lyout) %>%
  .$data %>% 
  mutate(spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
         support = as.numeric(label),
         support_class = cut(support, c(0,50,70,90,100)) %>% 
           as.character() %>% factor(levels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]"))
           )

p_tree <- (open_tree(
  ggtree(raxml_data, layout = lyout,
         aes(color = ifelse(clade == 0,
                            lab2spec(label),
                            clade2spec[as.character(clade)])), size = .25) %>%
    ggtree::rotate(200), 180))  +
  # geom_tippoint(size = .2) + 
  geom_tiplab2(aes(color = lab2spec(label), 
                   label = str_sub(
                     label, -6, -1)),
  size = GenomicOriginsScripts::plot_text_size_small / ggplot2:::.pt  *.6,#2.5, 
  hjust = -.1)+
  ggtree::geom_treescale(width = .002,
                         linesize = .2,
                         x = -.0007, y = 155, 
                         offset = -4,
                         fontsize = GenomicOriginsScripts::plot_text_size_small / ggplot2:::.pt,
                         color = clr_neutral) +
  xlim(c(-.0007,.0092)) +
  ggtree::geom_nodepoint(aes(fill = support_class, 
                             size = support_class),
                 shape = 21#, linewidth = 3
                 ) +
  scale_color_manual(values = c(ungrouped = clr_neutral, 
                                GenomicOriginsScripts::clr2),
                     guide = FALSE) +
  scale_fill_manual(values = c(`(0,50]` = "transparent",
                               `(50,70]` = "white",
                               `(70,90]` = "gray",
                               `(90,100]` = "black"),
                    drop = FALSE) +
  scale_size_manual(values = c(`(0,50]` = 0,
                               `(50,70]` = .4,
                               `(70,90]` = .4,
                               `(90,100]` = .4),
                    na.value = 0,
                    drop = FALSE)+
  guides(fill = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2,keyheight = unit(9,"pt")),
         size = guide_legend(title = "Node Support Class", title.position = "top", ncol = 2,keyheight = unit(9,"pt"))) +
  theme_void(base_size = GenomicOriginsScripts::plot_text_size_small  ) 

y_sep <- .05
x_shift <- -.03
p1 <- ggplot() +
  coord_equal(xlim = c(0, .93),
              ylim = c(-.01, .54),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(p_tree + theme(legend.position = "none")),
                    ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
                    xmin = -.1, xmax = 1.1) +
  annotation_custom(grob = cowplot::get_legend(p_tree),
                    ymin = .35, ymax = .54,
                    xmin = 0, xmax = .2) +
  theme_void()

data_ibd <- read_tsv(ibd_file) %>% 
  mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2)) 

set.seed(42)
p2 <- data_ibd %>% 
  as_tbl_graph() %>%
  mutate(spec = str_sub(name,-6,-4),
         loc = str_sub(name,-3,-1))  %>% 
  ggraph( layout = 'fr', weights = ibd_total) +
  geom_edge_link(aes(alpha = ibd_total), color = rgb(.1,.1,.1), edge_width = .15) +
  geom_node_point(aes(fill = spec,
                      shape = loc, color = after_scale(clr_darken(fill,.3))), size = 1.2) +
  scale_fill_manual("Species", values = GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))],
                    labels = GenomicOriginsScripts::sp_labs)+
  scale_edge_alpha_continuous(#limits = c(0,.1),
    range = c(0,1), guide = "none") +
  scale_shape_manual("Site", values = 21:23, labels = GenomicOriginsScripts::loc_names) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 21, size = 2.5)),
         shape = guide_legend(nrow = 2)) +
  coord_equal()  +
  theme(text = element_text(size = GenomicOriginsScripts::plot_text_size),
        panel.background = element_blank())

p_done <- (p1 + p2) / guide_area() +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, .07),
              guides = "collect") &
  theme(text = element_text(size = GenomicOriginsScripts::plot_text_size),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0)

hypo_save(plot = p_done,
          filename = "figures/SF9.png",
          width = GenomicOriginsScripts::f_width,
          height = GenomicOriginsScripts::f_width * .42,
          bg = "transparent",
          type = "cairo",
          dpi = 600,
          comment = plot_comment)

system("convert figures/SF9.png figures/SF9.pdf")
# hypo_save(plot = p_done,
#           filename = "figures/SF9.pdf",
#           width = 7.5 * scl,
#           height = 4 * scl,
#           device = cairo_pdf,
#           bg = "transparent",
#           comment = plot_comment)
