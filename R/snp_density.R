library(tidyverse)
library(hypogen)
library(ggpointdensity)
library(ggfx)
library(patchwork)

data <- vroom::vroom("2_analysis/revPoMo/window_stats.tsv.gz", delim = "\t") %>%
  left_join(hypo_chrom_start) %>%
  mutate(GPOS = GSTART + (START + END) / 2 )


p1 <- data %>%
  ggplot( aes( x = GPOS, y = SNP_density ) ) +
  geom_hypo_LG() +
  with_raster(geom_pointdensity(aes(group = CHROM),
                    size = .3)) +
  geom_smooth(aes(group = CHROM), color = "black", se = FALSE)

p2 <- data %>%
  ggplot( aes( x = GPOS, y = REL_COV ) ) +
  geom_hypo_LG() +
  with_raster(geom_pointdensity(aes(group = CHROM),
                                size = .3))+
  geom_smooth(aes(group = CHROM), color = "black", se = FALSE)

p1 / p2 &
  scale_x_hypo_LG() &
  scale_fill_hypo_LG_bg() &
  scale_color_viridis_c(option = "C") &
  theme_hypo() &
  theme(legend.position = "none")

ggsave("~/Desktop/snp_density.pdf", width = 16, height = 9, device = cairo_pdf)
