#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF12.R \
#     ressources/species_order_alpha.txt \
#     2_analysis/dstats/hyp_ld05_dtrios_BBAA.txt \
#     2_analysis/dstats/BBAA_ld05.csv \
#     2_analysis/dstats/BBAA_sign_ld05.csv
# ===============================================================
# This script produces Suppl. Figure 12 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("ressources/species_order_alpha.txt",
#           "2_analysis/dstats/hyp_ld05_dtrios_BBAA.txt",
#           "2_analysis/dstats/BBAA_ld05.csv",
#           "2_analysis/dstats/BBAA_sign_ld05.csv")
# script_name <- "R/fig/plot_SF12.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(prismatic)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
# config -----------------------
species_order_file <- as.character(args[1])
trios_file <- as.character(args[2])
bbaa_file <- as.character(args[3])
signif_file <- as.character(args[4])

two_chr_to_sorted_pair <- function(P2, P3, ...){
  n1 <- as.numeric(factor(P2, levels = sorter$group))
  n2 <- as.numeric(factor(P3, levels = sorter$group))
  char_sorted <- if(n1 < n2){c(P2, P3)} else { c(P3, P2)}
  str_c(char_sorted[[1]], "-", char_sorted[[2]])
}

sorter <- read_tsv(species_order_file, col_names = "group")

p_cap <- 16

data <- read_tsv(bbaa_file) %>%
  filter(p_adjusted <= .05) %>%
  mutate(pair = pmap_chr(.,two_chr_to_sorted_pair)) %>%
  group_by(pair) %>%
  mutate(p_is_min = p_adjusted == min(p_adjusted)) %>%
  filter(p_is_min) %>%
  mutate(d_is_max = Dstatistic == max(Dstatistic)) %>%
  filter(d_is_max) %>%
  ungroup() %>%
  separate(pair, into = c("p_left", "p_right"), sep = "-", remove = FALSE) %>%
  mutate(p_adjusted = if_else(p_adjusted < 10^-p_cap,10^-p_cap, p_adjusted ))

data_prep <- data %>%
  dplyr::select(p_left, p_right, Dstatistic, p_adjusted) %>%
  bind_rows(data %>%
              dplyr::select(p_left = p_right, p_right = p_left,
                            Dstatistic, p_adjusted))

data_full <- cross_df(list(p_left = sorter$group ,
                           p_right = sorter$group)) %>%
  left_join(data_prep) %>%
  mutate(p_left = factor(p_left, levels = sorter$group),
         p_right = factor(p_right, levels = rev(sorter$group)))

data_signif <- read_tsv(signif_file) %>%
  mutate(pair = pmap_chr(.,two_chr_to_sorted_pair)) %>%
  separate(pair, into = c("p_left", "p_right"), sep = "-", remove = FALSE) %>%
  mutate(p_left = factor(p_left, levels = sorter$group),
         p_right = factor(p_right, levels = rev(sorter$group))) %>%
  mutate(p_adjusted = if_else(p_adjusted < 10^-p_cap, 10^-p_cap, p_adjusted) )

clr <- scales::colour_ramp(colors = RColorBrewer::brewer.pal(5,"RdYlBu")[c(1,2,4,5)])((1:7)/7) %>% clr_saturate(.1)

d_lim <- c(0, .01)
p_lim <- c(1, p_cap)

p_done <- data_full %>%
  filter(p_left != "Outgroup" ) %>%
  mutate(check1 = as.numeric(p_left),
         check2 = as.numeric(p_right)) %>%
  filter(check2 < 17 - check1) %>%
  ggplot()+
  geom_tile(aes(x = p_left,
                y = p_right,
                fill = Dstatistic,
                color = after_scale(clr_darken(fill))) ) +
  geom_point(data = data_signif %>% filter(p_value < .05),
             aes(x = p_left,
                 y = p_right,
                 size = -log10(p_adjusted)),
             shape = 1,
             alpha = .4) +
  scale_fill_gradientn(colours = rev(clr),
                       limits = d_lim,
                       na.value = rgb(1,1,1,.2)) +
  scale_size(range = c(.1, 6),
             limits = c(1,16),
             breaks = c(1,8,16)) +
  guides(fill = guide_colorbar(title = "D",
                               title.position = "top",
                               barwidth = unit(.6,"npc"),
                               barheight = unit(4,"pt"),
                               order = 1),
         size = guide_legend(title = "-log<sub>10</sub> *( p<sub>adjusted</sub> )*",
                             order = 2,
                             title.position = "top")) +
  coord_equal() +
  theme_minimal(base_size = plot_text_size) +
  theme(legend.position = c(.95,1),
        legend.text.align = .5,legend.title.align = 1,
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.box.just = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5),
        axis.title = element_blank(),
        legend.title = element_markdown())

hypoimg::hypo_save(filename = "figures/SF12.pdf",
       width = f_width_half,
       height = f_width_half,
       device = cairo_pdf,
       bg = "transparent",
       comment = plot_comment)
