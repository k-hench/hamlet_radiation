#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF11.R 2_analysis/newhyb/nh_input/NH.Results/
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/newhyb/nh_input/NH.Results/")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(prismatic)
library(paletteer)
library(patchwork)
library(ggtext)

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
base_dir <- as.character(args[1])

folders <- dir(base_dir)

p_loc <- c("bel", "hon", "pan") %>%
  map(plot_loc)

# p_all <- p_loc[[1]] + 
#   p_loc[[2]] + 
#   p_loc[[3]] +
#   guide_area() +
#   plot_layout(ncol = 1,
#               heights = c(3,6,1.5,.5),
#               guides = "collect")

theme_hyb <-  function(legend.position = "none",...){
  list(scale_y_continuous(breaks = c(0,.5,1)),
       theme(legend.position = legend.position, 
        legend.background = element_rect(fill = "white",colour = rgb(1,1,1,0)),
        legend.direction = "horizontal",
        legend.justification = c(1,1),
        strip.text.y = element_markdown(angle = 0,hjust = 0), 
        ...))
}

label_spacer <- function(x, plus = 1.1){x + plus}


p <- (p_loc[[1]] +  guides(fill = guide_legend(title = "Hybrid Class")) + theme_hyb(legend.position = c(1,1)) ) + 
  (p_loc[[2]] + theme_hyb() ) + 
  (p_loc[[3]] + theme_hyb() )  + 
  plot_layout(ncol = 1, heights = c(10,15,3) %>% label_spacer())+ 
  plot_annotation(tag_levels = 'a')
  

ggsave(filename = "figures/SF11.pdf",
       plot = p,
       height = 16,
       width = 10, 
       device = cairo_pdf)

# get_data <- function (loc) 
# {
#   data <- map_dfr(.x = folders[str_detect(folders, loc)], .f = getPofZ, 
#                   base_dir = base_dir)
#   data 
# }
# 
# 
# 
# data_loc <- c("bel", "hon", "pan") %>%
#   map_dfr(get_data)
# 
# data_loc %>%
#   filter(!grepl(pattern = "_pure", bin))%>%
#   filter(post_prob > .99) %>%
#   arrange(IndivName) %>%
#   group_by(IndivName) %>%
#   nest()
