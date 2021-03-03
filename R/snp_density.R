library(tidyverse)
library(hypogen)
library(ggpointdensity)
library(ggfx)
library(patchwork)
library(prismatic)
library(ggstance)

data <- vroom::vroom("2_analysis/revPoMo/window_stats.tsv.gz", delim = "\t") %>%
  left_join(hypo_chrom_start) %>%
  mutate(GPOS = GSTART + (START + END) / 2 )

cov_tres <- quantile(data$REL_COV, probs = .66)
snp_tres <- quantile(data$SNP_density, probs = .66)

set.seed(42)
random_subset <- data %>%
  # filter(REL_COV > cov_tres,
  #        SNP_density > snp_tres) %>%
  sample_n(1000) 

p1 <- data %>%
  ggplot( aes( x = GPOS, y = SNP_density ) ) +
  geom_hypo_LG() +
  with_raster(geom_pointdensity(aes(group = CHROM),
                    size = .3, adjust = .1)) +
  geom_point(data = random_subset, shape = 1, color = "white") +
  geom_smooth(aes(group = CHROM), color = "black", se = FALSE) +
  scale_y_continuous(limits = c(0, .25)) +
  scale_color_viridis_c(option = "C", guide = FALSE) +
  scale_fill_hypo_LG_bg() +
  scale_x_hypo_LG() +
  theme_hypo()

p2 <- data %>%
  ggplot( aes( x = GPOS, y = REL_COV) ) +
  geom_hypo_LG() +
  with_raster(geom_pointdensity(aes(group = CHROM),
                                size = .3, adjust = .1)) +
  geom_point(data = random_subset, shape = 1, color = "white") +
  geom_smooth(aes(group = CHROM), color = "black", se = FALSE) +
  scale_color_viridis_c(option = "C",
                        guide = FALSE)+
  scale_fill_hypo_LG_bg() +
  scale_x_hypo_LG() +
  theme_hypo()

fll_vrds <- viridis::plasma(9) %>%
  clr_desaturate(shift = .3)

p3 <- data %>%
  ggplot() +
  geom_histogramh(aes(y = SNP_density,
                      x = ..count.. / 1000,
                      fill = ..count.., 
                      color = after_scale(clr_darken(fill))),
                  bins = 37) +
  geom_histogramh(data = random_subset, 
                  aes(x = ..count.. / 5,
                      y = SNP_density),
                  color = rgb(1,1,1,.6),
                  fill = rgb(1,1,1,.4),
                  bins = 37) +
  scale_fill_gradientn(colours = fll_vrds, guide = FALSE) +
  scale_y_continuous(limits = c(0, .25)) +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

p4 <- data %>%
  ggplot() +
  geom_histogramh(aes(y = REL_COV,
                      x = ..count.. / 1000,
                      fill = ..count.., 
                      color = after_scale(clr_darken(fill))),
                  bins = 37) +
  geom_histogramh(data = random_subset, 
                  aes(x = ..count.. / 5,
                      y = REL_COV),
                  color = rgb(1,1,1,.6),
                  fill = rgb(1,1,1,.4),
                  bins = 37) +
  scale_fill_gradientn(colours = fll_vrds, guide = FALSE) +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

p_done <-  p1 + p3 + p2 + p4 +
  plot_layout(widths = c(1, .07)) 

ggsave(plot = p_done,
       filename = "~/Desktop/snp_density.pdf",
       width = 16,
       height = 9,
       device = cairo_pdf)

random_subset %>%
  select(CHROM:END) %>%
  write_tsv(file = "2_analysis/revPoMo/random_1k_windows.bed.gz")