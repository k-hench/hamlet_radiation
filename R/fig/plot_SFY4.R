library(GenomicOriginsScripts)
library(patchwork)

clr_alt <- clr
clr_alt[["uni"]] <- "lightgray"
pca_fish_scale <- 1.15

plot_fish_lwd <- function (short, x = 0, y = 3, height = 5, width = 5, lwd = .15, line_color = "transparent") { 
  hypo_anno_l_lwd(sp_names[short], xmin = x - 0.5 * width, xmax = x + 
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 * 
                    height, lwd = lwd, line_color = line_color)
}
hypo_anno_l_lwd <- function (species, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, lwd = .15, line_color = "black") {
  stopifnot(length(species) == 1)
  stopifnot(is.character(species))
  stopifnot(species %in% hypo_img$spec)
  nr_species <- which(hypo_img$spec == species)
  annotation_custom(editGrob(grob = hypo_img$l[[nr_species]],
                             gPath = "GRID.picComplexPath.*", grep = TRUE,
                             gp = gpar( lwd = lwd, col = line_color
                             ), 
                             global = TRUE, strict = FALSE) ,
                    xmin = xmin, 
                    xmax = xmax, ymin = ymin, ymax = ymax)
}


pca_plot <- function(loc, mode, pc1 = 1, pc2 = 2){
  if(loc == ""){ loc2 <- "all" } else { loc2 <- loc }
  titles <- c(loc_names,"all") %>% set_names(nm = c("bel.", "hon.", "pan.", "flo", "all"))
  clr_pca <- c(clr_loc, "black") %>% set_names(nm = c("bel.", "hon.", "pan.", "flo", "all"))
  
  set.seed(42)
  evs <- str_c("2_analysis/pca/", loc , mode, ".exp_var.txt.gz") %>% 
    read_tsv()
  str_c("2_analysis/pca/", loc , mode,".scores.txt.gz") %>% 
    read_tsv() %>% 
    mutate(spec = str_sub(id, -6,-4)) %>%
    ggplot(aes_string(x = str_c("EV0",pc1), y = str_c("EV0",pc2), fill = "spec"))+
    ggforce::geom_mark_ellipse(aes(color = spec),
                               fill = "transparent",
                               linetype = 3,
                               size = .3,
                               expand = unit(5, "pt"))+
    geom_point(shape = 21, aes(color = after_scale(prismatic::clr_darken(fill))), size = .7) +
    # (pca_fish_pos$data[[ which(pca_fish_pos$loc == loc) ]] %>%
    #    pmap(plot_fish_lwd))+
    labs(x = str_c("PC",pc1," (", sprintf("%.1f",evs$exp_var[[ pc1 ]]), " %)"),
         y = str_c("PC",pc2," (", sprintf("%.1f",evs$exp_var[[ pc2 ]]), " %)"))+
    scale_fill_manual(values = clr)+
    scale_color_manual(values = clr_alt %>%
                         prismatic::clr_alpha(alpha = .7) %>%
                         set_names(nm = names(clr_alt)))+
    labs(title = str_c(titles[[ loc2 ]], " (", mode, ")"))+
    theme_minimal()+
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.background = element_rect(color = clr_pca[[ loc2 ]],
                                          fill = "transparent"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = plot_text_size),
          plot.title = element_text(color = clr_pca[[ loc2 ]]))
}

fish_tib <- tibble(short = names(clr),#[!names(clr), %in% c("flo", "tab", "tor")],
                   # x = ((seq_along(short) - 1)  %% 6 + 1)* 3 - 2.5,
                   # y = ((seq_along(short) - 1)  %/% 6 ) * -1.5 + .25,
                   x = (seq_along(short) - 1) * 3 + .5,
                   # x = c(0.5,  3.5,  7,  9.7, 12.25, 15.25, 18, 21.5
                         )

key_sz <- .75
sp_fam <- rep(c("H", "S", "H"), c(8, 2, 1)) %>%  set_names(nm = names(sp_names))
p_leg <- fish_tib %>% 
  ggplot() +
  # coord_equal(xlim = c(-.05, 12), expand = 0) +
  coord_equal(xlim = c(-.05, 33), 
              expand = 0) +
  geom_tile(aes(x = x, y = 0,
                fill = short, 
                color = after_scale(prismatic::clr_darken(fill, .25))),
            width = key_sz, height = key_sz, size = .3) +
  geom_text(aes(x = x + .6, y = 0,
                label = str_c(sp_fam[short], ". ", sp_names[short])), 
            hjust = 0, fontface = "italic", size = plot_text_size / ggplot2:::.pt) +
  # labs(subtitle = "Hamlet Species") +
  pmap(fish_tib, plot_fish_lwd, width = 1, height = 1, y = 0) +
  scale_fill_manual(values = clr, guide = FALSE) +
  theme_void()



p_done <- cowplot::plot_grid((tibble(loc = rep(c("bel.", "hon.", "pan.", ""), 4), 
                             mode = rep(c( "whg", "subset_non_diverged","whg", "subset_non_diverged"), each = 4),
                             pc1 = rep(c(1, 3), each = 8), 
                             pc2 = rep(c(2, 4), each = 8)) %>% 
                        pmap(pca_plot) %>% 
                        wrap_plots(ncol = 4) +
                        plot_annotation(tag_levels = "a") & 
                        theme(plot.background = element_blank())),
                     p_leg,
  ncol = 1 ,rel_heights = c(1,.05))

scl <- 1.3
hypo_save(p_done, filename = 'figures/SFY4.pdf',
          width = 9 * scl,
          height = 8.5 * scl,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
