### Tip-specific speciation rates across Fish Tree of Life
### ------------------------------------------------------

### Data from Rabosky et al. 2018 (Nature) -- Dryad: https://doi. org/10.5061/dryad.fc71cp4)
### Aug 11, 2021, Martin Helmkampf

library(GenomicOriginsScripts)

## Import FToL data
rates <- read_csv(file = "ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv")

p1 <- ggplot(rates, aes(lambda.tv)) +
  geom_histogram(binwidth = 0.25, colour="grey20", fill="grey80", size = .2) +
  labs(x = "Mean speciation rate", y = "Number of species") +
  geom_segment(aes(x = 2.5 , y = 200, xend = 2.5, yend = 1500), size = 0.2, color = "#0976BA") +
  annotate("text", x = 2.5, y = 1750, label = "Hamlets", size =  plot_text_size / ggplot2:::.pt, color = "#0976BA") +
  geom_segment(aes(x = 3.5 , y = 200, xend = 3.5, yend = 1500), size = 0.2, color = "grey60") +
  annotate("text", x = 3.5, y = 2150, label = "Haplo-\nchromines", size =  plot_text_size / ggplot2:::.pt, color = "grey60") +
  geom_segment(aes(x = 4.5 , y = 200, xend = 4.5, yend = 1500), size = 0.2, color = "grey60") +
  annotate("text", x = 4.5, y = 1750, label = "Labeobarbus", size = plot_text_size / ggplot2:::.pt, color = "grey60") +
  coord_cartesian(xlim = c(-.2, 5.1),
                  ylim = c(-200, 8400),
                  expand = 0) +
  theme_bw(base_size = plot_text_size) +
  theme(
    panel.grid.minor = element_blank())

ggsave("figures/SFx6.pdf", 
       plot = p1, device = cairo_pdf,
       width = f_width_half, height = f_width_half*.6)

