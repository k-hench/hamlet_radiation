#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFX.R
# ===============================================================
# This script produces Suppl. Figure X of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/summaries/fst_outliers_998.tsv", "2_analysis/geva/", "2_analysis/GxP/bySNP/" )
# script_name <- "R/fig/plot_SFX.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
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

import_geva_data <- function(chrom, start, end, ...){
  vroom::vroom(file = str_c(geva_path, chrom,".sites.txt.gz"), delim = " ") %>%
    left_join(vroom::vroom(str_c(geva_path, chrom,".marker.txt.gz"), delim = " ")) %>%
    arrange(Position) %>%
    filter(between(Position,start,end)) %>%
    select(Chromosome,Position,MarkerID, Clock, Filtered:PostMedian) %>%
    mutate(CHROM = str_c("LG", str_pad(Chromosome,width = 2,pad = "0")))
}

gxp_importer <- function(trait, chrom, start, end ){
  vroom::vroom(str_c(gxp_path,trait,".lm.GxP.txt.gz"), delim = "\t") %>%
    filter(CHROM == chrom, between(POS, start, end) ) %>%
    select(CHROM, POS, p_wald) %>%
    set_names(nm = c("CHROM", "Position", trait))
}

import_gxp_data <- function(chrom, start, end, ...){
  bar_data <- gxp_importer("Bars", chrom, start, end)
  peduncle_data <- gxp_importer("Peduncle", chrom, start, end)
  snout_data <-gxp_importer("Snout", chrom, start, end)
  
  bar_data %>%
    left_join(peduncle_data, by = c(CHROM = "CHROM", Position = "Position")) %>%
    left_join(snout_data, by = c(CHROM = "CHROM", Position = "Position"))
}

get_gxp_and_geva <- function(gid, chrom, start, end, ...){
  import_geva_data(chrom, start, end) %>%
    left_join(import_gxp_data(chrom, start, end)) %>%
    mutate(gid = gid)
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

grid_piece_x <- function(g, x1, x2, ...){
  annotation_custom(grob = g, xmin = x1, xmax = x2, ...)
}

theme_gradient_bg_x <- function(xlim_in, colors = c("#7B0664", "#E32219"), fun = identity, step = 1,...){
  g <- rasterGrob(cbind(colors[[1]], colors[[2]]), width = unit(1, "npc"),
                  height = unit(1, "npc"), 
                  interpolate = TRUE)
  
  xstart <- floor(fun(min(xlim_in)))
  xend <- ceiling(fun(max(xlim_in)))
  
  breaks <- seq(xstart,xend, by = step)
  
  purrr::map2(.x = breaks[1:(length(breaks)-1)],
              .y = breaks[2:length(breaks)],
              .f = grid_piece_x, 
              g = g, ...)
}

data <- outlier_data[c(2,13,14),] %>%
  pmap_dfr(get_gxp_and_geva)

xrange <- c(100,10^6)
color <- rgb(1, 0.5, 0.16)

base_length <- 8
base_lwd <- .15
base_line_clr <- "black"

splitage <- tibble(intercept = 5000)

p <- data %>%
  pivot_longer(names_to = "trait",
               values_to = "p_wald",
               cols = Bars:Snout) %>%
  filter(Clock == "J", Filtered == 1) %>%
  ggplot(aes(x = PostMedian,y = p_wald)) +
  #geom_vline(data = splitage, aes(xintercept = intercept), linetype = 3, color = rgb(0,0,0,.9)) +
  theme_gradient_bg_x(xlim_in = xrange, #range(economics$unemploy),
                      colors = c( rgb(.95,.95,.95,0),
                                  rgb(.1,.1,.1,.15)),
                      fun = log10) +
  geom_pointdensity(size = .75)+
  facet_grid(gid ~ trait, scales = "free_y")+
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans=reverselog_trans(10),
                      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#  scale_color_viridis_c(option = "B")+
#  scico::scale_color_scico(palette = "berlin")+
  scale_color_gradientn(colours = scico::scico(n = 9,palette = "berlin") %>% clr_desaturate(shift = .1) #%>% clr_darken(shift = .25)
                        )+
  labs(y = "*p* value <sub>Wald</sub>",
       x  = "Allele Age (gen., Post Median)")+
  guides(color = guide_colorbar(barwidth = unit(150, "pt"),
                                barheight = unit(7, "pt")))+
  theme_minimal()+
  theme(#text = element_text(family = "Futura Lt BT"),
        axis.title.y = element_markdown(),
        legend.position = "bottom",
        plot.subtitle = element_markdown(),
        axis.line = element_line(colour = base_line_clr,
                                 size = base_lwd), 
        strip.background = element_rect(fill = rgb(.95,.95,.95),colour = base_line_clr,size = base_lwd),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

scl <- 1
hypo_save(plot = p,
          filename = "figures/SFX.pdf",
          width = 10 * scl,height = 7 * scl,
          comment = plot_comment,
          device = cairo_pdf,
          bg = "transparent")

