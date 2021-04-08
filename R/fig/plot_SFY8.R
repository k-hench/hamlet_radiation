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
           "lg12.4_155N.raxml.support",
           "lg04.1_hySN.raxml.support",
           "lg12.3_hySN.raxml.support",
           "lg12.4_hySN.raxml.support")

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
  mutate(outgroup = list("PL17_160floflo", "PL17_160floflo", "PL17_160floflo",
                         c("28393torpan", "s_tort_3torpan", "20478tabhon"),
                         c("28393torpan", "s_tort_3torpan", "20478tabhon"),
                         c("28393torpan", "s_tort_3torpan", "20478tabhon")),
         rooted_tree = map2(.x = tree,
                            .y =  outgroup,
                            .f = function(tree, outgroup){ape::root(phy = tree, outgroup = outgroup)}),
         tree_data = map(.x = rooted_tree, get_tree_data))

plot_tree_clasic <- function(tree, trait = "Bars", outlier = NULL, hjust = 1){
  tree_data <- get_tree_data(tree)
  ggtree(tree_data, mapping = aes(color = spec)) +
    # annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == trait][[1]],
    #                   xmin = 0, xmax = .001,
    #                   ymin = 145, ymax = 155)+
    geom_tiplab(aes(label = str_sub(label, -6, -1)),
                size = 2, hjust = hjust) +
    # geom_tippoint(size = .4) +
    geom_nodepoint(data = tree_data %>% 
                     filter(!isTip, !is.na(as.numeric(label))),
                   aes(fill =  as.numeric(label)),
                   size = .85,
                   shape = 21) + 
    ggtitle(outlier) +
    # scale_x_continuous(limits = c(0,.00485)) +
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


p1 <- plot_tree_clasic(data$tree_data[[1]], outlier = data$outlier[[1]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[4]], outlier = data$outlier[4], trait = "Snout") +
  plot_layout(guides = "collect") +
  plot_annotation(title = data$outlier[[1]])

p2 <- plot_tree_clasic(data$tree_data[[2]], outlier = data$outlier[[2]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[5]], outlier = data$outlier[5], trait = "Snout") +
  plot_layout(guides = "collect")+
  plot_annotation(title = data$outlier[[2]])

p3 <- plot_tree_clasic(data$tree_data[[3]], outlier = data$outlier[[3]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[6]], outlier = data$outlier[6], trait = "Snout") +
  plot_layout(guides = "collect")+
  plot_annotation(title = data$outlier[[3]])


p1 /
  p2 /
  p3 +
  plot_layout(guides = "collect")

p_done <- plot_tree_clasic(data$tree_data[[1]], outlier = data$outlier[[1]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[4]], outlier = "", trait = "Snout", hjust = -.1) +
  plot_tree_clasic(data$tree_data[[2]], outlier = data$outlier[[2]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[5]], outlier = "", trait = "Snout", hjust = -.1) +
  plot_tree_clasic(data$tree_data[[3]], outlier = data$outlier[[3]], trait = "Snout") +
  scale_x_reverse() +
  plot_tree_clasic(data$rooted_tree[[6]], outlier = "", trait = "Snout", hjust = -.1) +
  plot_layout(guides = "collect", ncol = 2, widths = c(1,.6))


scl <- 1.6
ggsave(plot = p_done,
       filename = "figures/SFY8.pdf",
       width = f_width * scl,
       height = f_width * scl * 1.2 ,
       device = cairo_pdf)


p_classic <- plot_tree_clasic(data$tree_data[[1]], outlier = data$outlier[[1]], trait = "Snout", hjust = -.1) +
  plot_tree_clasic(data$tree_data[[2]], outlier = data$outlier[[2]], trait = "Bars", hjust = -.1) +
  plot_tree_clasic(data$tree_data[[3]], outlier = data$outlier[[3]], trait = "Peduncle", hjust = -.1) +
  plot_tree_clasic(data$rooted_tree[[4]], outlier = data$outlier[4], trait = "Snout", hjust = -.1) +
  plot_tree_clasic(data$rooted_tree[[5]], outlier = data$outlier[[5]], trait = "Bars", hjust = -.1) +
  plot_tree_clasic(data$rooted_tree[[6]], outlier = data$outlier[[6]], trait = "Peduncle", hjust = -.1) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", 
              nrow = 2,
              widths = c(1,1,1))

scl <- 1.6
ggsave(plot = p_classic,
       filename = "figures/SFY8b.pdf",
       width = f_width * scl,
       height = f_width * .7 *scl,
       device = cairo_pdf)
