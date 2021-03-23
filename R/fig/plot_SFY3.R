#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFY3.R ~/work/puebla_lab/stash/hyp155_n_0.33_mac4_5kb.raxml.support.bs-tbe
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("~/work/puebla_lab/stash/hyp155_n_0.33_mac4_5kb.raxml.support.bs-tbe")
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(ape)
library(ggtree)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
tree_file <- as.character(args[1])

raxml_tree <- read.tree(tree_file) 
raxml_tree_rooted <- root(phy = raxml_tree, outgroup="PL17_160floflo")
clr_neutral <- rgb(.6, .6, .6)
lyout <- 'circular'

lab2spec <- function(label){
  x <- str_sub(label, start = -6, end = -4) %>% str_remove(.,"[0-9.]{1,3}$") %>% str_remove(.," ")
  ifelse(x == "",'ungrouped', x)
}

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
  ggtree::rotate(200) %>%
  .$data 

p_tree <- (open_tree(
  ggtree(raxml_tree_rooted_grouped, layout = lyout,
         aes(color = ifelse(clade == 0,
                            lab2spec(label),
                            clade2spec[as.character(clade)]))) %>%
    ggtree::rotate(200), 180))  +
  # geom_tree(layout = lyout, aes(color = lab2spec(label)), size = .1) +
  geom_tippoint(size = .4) + 
  geom_tiplab2(aes(color = lab2spec(label), 
                   label = str_sub(
                     label, -6, -1)
  ),
  size = 3, hjust = -.1)+
  # add node-support layer
  # geom_nodelab(data = raxml_data %>% filter(!isTip, !is.na(as.numeric(label))),
  #                        aes(label = sprintf("%.2f", as.numeric(label))),
  #              size = 3,
  #              color = "black") +
  geom_nodepoint(data = raxml_data %>% filter(!isTip, !is.na(as.numeric(label))),
                 aes(fill =  as.numeric(label)),
                 size = 1.5,
                 shape = 21) +
  # geom_nodelab(data = raxml_data %>% filter(!isTip, !is.na(as.numeric(label))),
  #              aes(label = node),
  #              size = 3,
  #              color = "black") +
  ggtree::geom_treescale(width = .002,
                         x = -.0007, y = 155, 
                         offset = -3,fontsize = 3,
                         color = clr_neutral) +
  xlim(c(-.0007,.0092)) +
  scale_color_manual(values = c(ungrouped = clr_neutral, GenomicOriginsScripts::clr2),
                     guide = FALSE) +
  # scale_fill_distiller(palette = "Greys",
  #                      direction = 1,
  #                      limits = c(0, 1),
  #                      breaks = c(0, .5, 1),
  #                      labels = c("0","0.5", "1")) +
  # scale_fill_gradientn(colours = c(rev(RColorBrewer::brewer.pal(3,"Greens")[c(2,3)]),
  #                                  rev(RColorBrewer::brewer.pal(4, "Greys"))),
  #                      values = c(1, .7, seq(.7, 0,length.out = 4)),
  #                      limits = c(0, 1),
  #                      breaks = c(0, .5, .7, 1),
  #                      labels = c("0","0.5","0.7", "1")) +
scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Greys")[c(1,1,4)]),
                     values = c(1, .5, 0),
                     limits = c(0, 1),
                     breaks = c(0, .5, 1),
                     labels = c("0","0.5", "1")) +
  guides(fill = guide_colorbar(title = "Node Support",
                               title.position = "top",
                               barwidth = unit(87, "pt"),
                               barheight = unit(4, "pt"),
                               direction = "horizontal",
                               label.position = "top",
                               ticks.colour = "black")) +
  theme_void() 


p_done <- ggplot() +
  coord_equal(xlim = c(0, .93), 
              ylim = c(-.04, .52),
              expand = 0) +
  # annotation_custom(grob = ggplotGrob(p_tree + theme(legend.position = "none")),
  #                   ymin = -.575, ymax = .575,
  #                   xmin = -.075, xmax = 1.075) +
  annotation_custom(grob = ggplotGrob(p_tree + theme(legend.position = "none")),
                    ymin = -.6, ymax = .6,
                    xmin = -.1, xmax = 1.1) +
  annotation_custom(grob = cowplot::get_legend(p_tree), 
                    ymin = 0, ymax = .075,
                    xmin = 0, xmax = .907) +
  theme_void()

scl <- 1.5
hypo_save(plot = p_done,
       filename = "figures/SFY3.pdf", 
       width = 7.5 * scl,
       height = 4.5 * scl, 
       device = cairo_pdf, 
       bg = "transparent",
       comment = plot_comment)
