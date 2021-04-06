#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFY7.R ~/work/puebla_lab/stash/
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("~/work/puebla_lab/stash/")
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)
library(patchwork)
library(ggplotify)
library(phangorn)
cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
tree_path <- as.character(args[1])
trees <- c("lg04.1_155N.raxml.support",
           "lg12.3_155N.raxml.support",
           "lg12.4_155N.raxml.support")

import_tree <- function(file){
  tag = file %>% str_remove(".*/") %>%
    str_remove("N.raxml.support")
  tibble(outlier = str_sub(tag, 1, 6) %>% 
           str_to_upper() %>% 
           str_replace("\\.","_"),
         tree = list(read.tree(file)))
}

get_tree_data <- function(tree){
  ggtree(tree, layout = "circular") %>% 
    .$data %>% 
    mutate(spec = ifelse(isTip,  str_sub(label, -6, -4), "ungrouped"))
}
clr_neutral <- rgb(.6, .6, .6)
data <- str_c(tree_path, trees) %>% 
  map_dfr(import_tree) %>% 
  mutate(rooted_tree = map(.x = tree, .f = root, outgroup = "PL17_160floflo"),
         tree_data = map(.x = rooted_tree, get_tree_data))

plot_tree <- function(tree_data,xlim =  c(-.0005,.00475)){
  open_tree(ggtree(tree_data, layout = "circular",
                   aes(color = spec)), 180) +
    geom_tiplab(aes(label = str_sub(label, -6, -1)),
                size = 2, hjust = -.1) +
    geom_tippoint(size = .4) +
    geom_nodepoint(data = tree_data %>% 
                     filter(!isTip, !is.na(as.numeric(label))),
                   aes(fill =  as.numeric(label)),
                   size = 1.5,
                   shape = 21) +
    scale_color_manual(values = c(ungrouped = clr_neutral,
                                  GenomicOriginsScripts::clr2),
                       guide = FALSE) +
    scale_x_continuous(limits = xlim)+
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
}

p1 <- ggplot() +
  coord_equal(xlim = c(0, .93), 
              ylim = c(-.04, .52),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Snout"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = 0.02, ymax = .1)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(data$tree_data[[1]]) + theme(legend.position = "none")),
                    ymin = -.6, ymax = .7,
                    xmin = -.1, xmax = 1.2) +
  ggtitle(data$outlier[[1]]) +
  theme_void()+
  theme(plot.title = element_text(hjust = .5))

p2 <- ggplot() +
  coord_equal(xlim = c(0, .93), 
              ylim = c(-.04, .52),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Bars"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = 0.02, ymax = .1)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(data$tree_data[[2]],xlim = c(-.0005,.0045)) + 
      theme(legend.position = "none")),
    ymin = -.645, ymax = .76,
    xmin = -.17, xmax = 1.2) +
  ggtitle(data$outlier[[2]]) +
  theme_void()+
  theme(plot.title = element_text(hjust = .5))

p3 <- ggplot() +
  coord_equal(xlim = c(0, .93), 
              ylim = c(-.04, .52),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Peduncle"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = 0.02, ymax = .1)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(data$tree_data[[3]],xlim = c(-.0005,.0045)) + 
      theme(legend.position = "none")),
    ymin = -.645, ymax = .76,
    xmin = -.14, xmax = 1.22) +
  ggtitle(data$outlier[[3]]) +
  theme_void()+
  theme(plot.title = element_text(hjust = .5))

p_done <- plot_grid(p1, p2, p3,
          cowplot::get_legend( plot_tree(data$tree_data[[3]])),
          ncol = 1,
          rel_heights = c(1,1,1,.2), 
          labels = c(letters[1:3], ""),
          label_fontface = "plain")

scl <- 1
ggsave(plot = p_done,
       filename = "figures/SFY7.pdf",
       width = f_width * scl,
       height = f_width * 1.7 *scl,
       device = cairo_pdf)

plot_tree_clasic <- function(tree_data, trait = "Bars", outlier = NULL){
ggtree(tree_data, mapping = aes(color = spec)) +
    annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == trait][[1]],
                      xmin = 0, xmax = .001,
                      ymin = 145, ymax = 155)+
  geom_tiplab(aes(label = str_sub(label, -6, -1)),
              size = 2, hjust = -.1) +
  geom_tippoint(size = .4) +
  geom_nodepoint(data = tree_data %>% 
                   filter(!isTip, !is.na(as.numeric(label))),
                 aes(fill =  as.numeric(label)),
                 size = .85,
                 shape = 21) + 
  ggtitle(outlier) +
  scale_x_continuous(limits = c(0,.00485)) +
  scale_color_manual(values = c(ungrouped = clr_neutral,
                                GenomicOriginsScripts::clr2),
                     guide = FALSE) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Greys")[c(1,1,4)]),
                       values = c(1, .5, 0),
                       limits = c(0, 1),
                       breaks = c(0, .5, 1),
                       labels = c("0","0.5", "1")) +
  guides(fill = guide_colorbar(title = "Support",
                               title.position = "top",
                               barheight = unit(87, "pt"),
                               barwidth = unit(4, "pt"),
                               direction = "vertical",
                               ticks.colour = "black")) +
  theme_void() 
}

p_classic <- plot_tree_clasic(data$tree_data[[1]], outlier = data$outlier[[1]], trait = "Snout") +
  plot_tree_clasic(data$tree_data[[2]],outlier = data$outlier[[2]], trait = "Bars") +
  plot_tree_clasic(data$tree_data[[3]],outlier = data$outlier[[3]], trait = "Peduncle") +
  plot_annotation(tag_levels = "a") +
  guide_area() +
  plot_layout(guides = "collect", 
              nrow = 1,
              widths = c(1,1,1,.2))

scl <- 1.3
ggsave(plot = p_classic,
       filename = "figures/SFY7_classic.pdf",
       width = f_width * scl,
       height = f_width * 1 *scl,
       device = cairo_pdf)
