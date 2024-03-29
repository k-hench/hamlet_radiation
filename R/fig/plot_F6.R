#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F6.R \
    # 2_analysis/summaries/fst_outliers_998.tsv \
    # 2_analysis/geva/ \
    # 2_analysis/GxP/bySNP/
# ===============================================================
# This script produces Figure 6 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/summaries/fst_outliers_998.tsv",
#            "2_analysis/geva/", "2_analysis/GxP/bySNP/" )
# script_name <- "R/fig/plot_F6.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ggtext)
library(ggpointdensity)
library(scales)
library(grid)
library(prismatic)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name)))
args = args[7:length(args)]
cat(' ')
cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args),': ',crayon::green(args),'\n'))
cli::rule(right = getwd())

# config -----------------------
outlier_file <- as.character(args[1])
geva_path <- as.character(args[2])
gxp_path <- as.character(args[3])

outlier_data <- read_tsv(outlier_file)

data <- outlier_data[c(2, 13, 14),] %>%
  pmap_dfr(get_gxp_and_geva)

xrange <- c(100, 10^6)
color <- rgb(1, 0.5, 0.16)

base_length <- 8
base_lwd <- .15
base_line_clr <- "black"

splitage <- tibble(intercept = 5000)

gid_label <- c( LG04_1 = "LG04 (A)", LG12_3 = "LG12 (B)", LG12_4 = "LG12 (C)" )

gxp_clr <- c(Bars = "#79009f", Snout = "#E48A00", Peduncle = "#5B9E2D") %>%
  darken(factor = .95) %>%
  set_names(., nm = c("Bars", "Snout", "Peduncle"))

annotation_grobs <- tibble(svg = hypo_trait_img$grob_circle[hypo_trait_img$trait %in% c( 'Snout', 'Bars', 'Peduncle')],
                           layer = c(4,3,7),
                           color = gxp_clr[c(1,3,2)]) %>%
  purrr::pmap(.l = ., .f = hypo_recolor_svg) %>%
  set_names(nm = c( "LG12_3","LG12_4","LG04_1"))

annotation_grobs$LG12_3 <- hypo_recolor_svg(annotation_grobs$LG12_3,
                                            layer = 7, color = gxp_clr[[1]] %>% 
                                              clr_desaturate %>% clr_lighten(.25))

annotation_grobs_tib <- tibble(gid = names(annotation_grobs),
                               grob = annotation_grobs) %>%
  mutate( gid_label = gid_label[gid],
          trait = factor( c( "Bars", "Peduncle", "Snout"),
                          levels = c("Snout", "Bars", "Peduncle")))

highlight_rects <- tibble(trait = factor( c("Snout", "Bars", "Peduncle"),
                                          levels = c("Snout", "Bars", "Peduncle")),
                          gid_label = gid_label)

p_done <- data %>%
  pivot_longer(names_to = "trait",
               values_to = "p_wald",
               cols = Bars:Snout) %>%
  mutate(trait = factor(trait, levels = c("Snout", "Bars", "Peduncle")),
         gid_label = gid_label[gid]) %>%
  filter(Clock == "J",
         Filtered == 1) %>%
  ggplot() +
  geom_rect(data = highlight_rects,
            aes( xmin = 0, xmax = Inf,
                 ymin = 0, ymax = Inf),
            color = rgb(.75,.75,.75),
            size = .4,
            fill = rgb(.9,.9,.9,.5))+ 
  hypoimg::geom_hypo_grob(inherit.aes = FALSE,
                          data = annotation_grobs_tib,
                          aes(grob = grob), x = .15,  y = .78, angle = 0, width = .35, height =.35)+
  geom_pointdensity(size = plot_size,
                    aes(x = PostMedian,y = p_wald))+
  facet_grid(gid_label ~ trait, scales = "free_y")+
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans = reverselog_trans(10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_viridis_c("Density",  option = "B")+
  labs(y = "G x P *p* value <sub>Wald</sub>",
       x  = "Derived allele age (generations)")+
  guides(color = guide_colorbar(barwidth = unit(120, "pt"),
                                barheight = unit(3, "pt")))+
  theme_minimal()+
  theme(text = element_text(size = plot_text_size),
        axis.title.y = element_markdown(),
        legend.position = "bottom",
        plot.subtitle = element_markdown(),
        axis.line = element_line(colour = base_line_clr,
                                 size = base_lwd), 
        strip.background = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = plot_lwd)
  )

hypo_save(plot = p_done,
          filename = "figures/F6.pdf",
          width = f_width_half,
          height = f_width_half,
          comment = plot_comment,
          device = cairo_pdf,
          bg = "transparent")
