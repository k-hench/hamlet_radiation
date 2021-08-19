library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ape)
library(ggtree)
library(tidygraph)
library(ggraph)
library(patchwork)

# config -----------------------
tree_hypo_file <- "~/Downloads/hyp155_a_0.33_mac4_5kb.raxml.support"

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
         aes(color = lab2spec(label)), size = .25) %>%
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
              ylim = c(-.01, .4),
              expand = 0) +
  annotation_custom(grob = ggplotGrob(p_tree + theme(legend.position = "none")),
                    ymin = -.6 + (.5 * y_sep), ymax = .6 + (.5 * y_sep),
                    xmin = -.1, xmax = 1.1) +
  annotation_custom(grob = cowplot::get_legend(p_tree),
                    ymin = .3, ymax = .4,
                    xmin = 0, xmax = .2) +
  theme_void()
p1


ggsave(plot = p1,
          filename = "figures/SFx8.pdf",
          width = GenomicOriginsScripts::f_width,
          height = GenomicOriginsScripts::f_width * .5,
          bg = "transparent",
          device = cairo_pdf)
