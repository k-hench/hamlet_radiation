#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFY6.R ~/work/puebla_lab/stash/pomo/
# ===============================================================
# This script produces Figure 1 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("~/work/puebla_lab/stash/pomo/")
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
pomo_path <- as.character(args[1])
trees <- dir(pomo_path)

import_tree <- function(file){
  tag = file %>% str_remove(".*/") %>%
    str_remove("_pop.cf.treefile")
  tibble(outlier = str_sub(tag, 1, 6) %>% str_to_upper() %>% str_replace("\\.","_"),
         type = str_sub(tag, 8, -1) ,
         tree = list(read.tree(file)))
}

root_hamlets <- function(tree, type){
  outgr <- switch (type,
                   "155" = "floflo",
                   "hyS" = c("tabhon", "torpan")
  )
  root(phy = tree, outgroup = outgr)  
}

data <- str_c(pomo_path, trees) %>% 
  map_dfr(import_tree) %>% 
  mutate(rooted_tree = map2(.x = tree,.y = type, .f = root_hamlets))

conditional_rotate <- function(tree, rotate = TRUE, root_node){
  if(rotate){tree %>% ggtree::rotate(root_node)} else{tree}
}

conditional_highlight <- function(tree,
                                  highl = TRUE, higl_node, 
                                  support_guide = FALSE ){
  if(highl){
    p <- tree %>% 
      ggtree::groupClade(.node = higl_node, group_name =  "clade") %>% 
      ggtree(layout = "circular", aes(color = clade)) 
  } else {
      p <- tree %>% 
      ggtree(layout = "circular")
  }
  p <- p +
    geom_nodepoint(aes(fill =  support_class, size = support_class),
                   shape = 21) +
    scale_fill_manual(values = c(`(0,0.5]` = "transparent",
                                 `(0.5,0.7]` = "white",
                                 `(0.7,0.9]` = "gray",
                                 `(0.9,1]` = "black"),
                      drop = FALSE,
                      guide = FALSE) +
    scale_size_manual(values = c(`(0,0.5]` = 0,
                                 `(0.5,0.7]` = 1.5,
                                 `(0.7,0.9]` = 1.5,
                                 `(0.9,1]` = 1.5),
                      na.value = 0,
                      drop = FALSE,
                      guide = FALSE) 
  
  if(support_guide){p <- p +   guides(fill = guide_legend(title = "Support", title.position = "top", ncol = 1),
                                      size = guide_legend(title = "Support",title.position = "top", ncol = 1)) }
  p
  }

plot_tree <- function(tree, angle_in = 0,
                      color = "red", 
                      higl_clade = TRUE, higl_node = NA,
                      root_node = NA, rotate_root = FALSE, 
                      xlim = c(-2,4.5), 
                      support_guide = FALSE){
  (tree %>% 
     ggtree(layout = "circular") %>% 
     .$data %>% 
     mutate(support = as.numeric(label) /100,
            support_class = cut(support, c(0,.5,.7,.9,1)) %>% 
              as.character() %>%
              factor(levels = c("(0,0.5]", "(0.5,0.7]", "(0.7,0.9]", "(0.9,1]"))
     ) %>% 
     conditional_highlight(highl = higl_clade,higl_node = higl_node, support_guide = support_guide) %>% 
     conditional_rotate(rotate = rotate_root, root_node = root_node) %>% 
     open_tree(angle =  180) %>% 
     rotate_tree(angle = angle_in)) +
    geom_tiplab(offset = diff(xlim) * .07) +
    geom_tippoint(size = .4) +
    scale_color_manual(values  = c(`0` = "black", `1` = color),
                       guide = FALSE) +
    scale_x_continuous(limits = xlim) +
    theme_void()
  
}

# plot_tree(data$rooted_tree[[1]], higl_node = 28, color = "blue", xlim = c(-.05,.07)) +
#   geom_node_text(aes(label = node), color ="red") +
#   theme_gray()
# 
# plot_tree(data$rooted_tree[[2]], higl_node = 28, color = "blue") +
#   geom_node_text(aes(label = node), color ="red") +
#   theme_gray()
# 
# 

twisst_clr <- c(Blue = "#0140E5", Bars = "#E32210", Butter = "#E4E42E")

yshift <-  0#.025
p1 <- ggplot() +
  coord_equal(xlim = c(0, 1), 
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Snout"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = .4, ymax = .6)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(midpoint(data$rooted_tree[[2]]),
                                                higl_node = 24,# 32,
                                                color = twisst_clr[["Blue"]],
                                                root_node = 18,
                                                rotate_root = TRUE,
                                                angle_in = 180,
                                                xlim = c(-.8, 2))
    ),
                    ymin = 0 - yshift, ymax = 1 - yshift,
                    xmin = 0, xmax = 1) +
  annotation_custom(grob = ggplotGrob(
    plot_tree(midpoint(data$tree[[1]]), 
              higl_node = 21,
              angle_in = -3,
              color = twisst_clr[["Blue"]], 
              xlim = c(-.045, .11))
    ),
                    ymin = 0 + yshift, ymax = 1 + yshift,
                    xmin = 0, xmax = 1) +
  ggtitle(data$outlier[[1]]) +
  theme_void() +
  theme(plot.title = element_text(hjust = .5))

p2 <- ggplot() +
  coord_equal(xlim = c(0, 1), 
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Bars"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = .4, ymax = .6)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(
      midpoint(data$rooted_tree[[4]]),
      higl_node = 31,
      color = twisst_clr[["Bars"]],
      root_node = 18,
      rotate_root = TRUE,
      angle_in = 180,
      xlim = c(-.54, 1.3)) 
    ),
                    ymin = 0 - yshift, ymax = 1 - yshift,
                    xmin = 0, xmax = 1) +
  annotation_custom(grob = ggplotGrob(
    plot_tree(midpoint(data$tree[[3]]), 
              higl_node = 27,
              color = twisst_clr[["Bars"]], 
              xlim = c(-.035,.08))
  ),
                    ymin = 0 + yshift, ymax = 1 + yshift,
                    xmin = 0, xmax = 1) +
  ggtitle(data$outlier[[3]]) +
  theme_void() +
  theme(plot.title = element_text(hjust = .5))

p3 <- ggplot() +
  coord_equal(xlim = c(0, 1), 
              ylim = c(0, 1),
              expand = 0) +
  annotation_custom(grob = hypo_trait_img$grob_circle[hypo_trait_img$trait == "Peduncle"][[1]],
                    xmin = .4, xmax = .6,
                    ymin = .4, ymax = .6)+
  annotation_custom(grob = ggplotGrob(
    plot_tree(midpoint(data$rooted_tree[[6]]),
                                                higl_node = 32,#26,
                                                color = twisst_clr[["Butter"]],
                                                root_node = 18,
                                                rotate_root = TRUE,
                                                angle_in = 180,
                                                xlim = c(-.7, 1.7))
                                      ),
                    ymin = 0 - yshift, ymax = 1 - yshift,
                    xmin = 0, xmax = 1) +
  annotation_custom(grob = ggplotGrob(
    plot_tree(midpoint(data$tree[[5]]), 
              higl_node = 28,
              color = twisst_clr[["Butter"]], 
              xlim = c(-.03, .07))
    ),
                    ymin = 0 + yshift, ymax = 1 + yshift,
                    xmin = 0, xmax = 1) +
  ggtitle(data$outlier[[5]]) +
  theme_void() +
  theme(plot.title = element_text(hjust = .5))

p_done <- plot_grid(p1, p2, p3,
                    get_legend(plot_tree(data$rooted_tree[[5]], 
                                         higl_node = 23,
                                         color = twisst_clr[["Butter"]], 
                                         xlim = c(-.01,.136),
                                         support_guide = TRUE)),
                    nrow = 1,
                    rel_widths = c(1,1,1,.2), labels = c(letters[1:3], ""),
                    label_fontface = "plain")
scl <- 2.1
ggsave(plot = p_done,
       filename = "figures/SFY6.pdf",
       width = f_width * scl,
       height = f_width * .36 * scl,
       device = cairo_pdf)
