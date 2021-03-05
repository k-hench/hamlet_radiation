library(GenomicOriginsScripts)
library(prismatic)

sort_element <- function(x1, x2, idx){
 sort(c(x1, x2))[[idx]] 
}

data <- read_tsv("2_analysis/fst_signif/genepop_summary.tsv.gz") %>%
  separate(Population, into = c("pop1", "pop2"), sep = "-") %>%
  mutate(p = str_replace(df, pattern = "<", replacement =  ""),
         p = str_replace(p, pattern = "Highly_sign.",
                         replacement =  "0") %>%
           as.numeric(),
         p1 = map2_chr(pop1,pop2,sort_element,idx = 1),
         p2 = map2_chr(pop1,pop2,sort_element,idx = 2))

bg_black <- rgb(0, 0, 0, .75)

data %>%
  ggplot(aes(x = p1, y = p2, fill = p))+
  coord_equal() +
  ggfx::with_blur(geom_tile(color = bg_black,
                            fill = bg_black,
            width = .8, height = .8),
            sigma = 1)+
  geom_tile(aes(color = after_scale(prismatic::clr_darken(fill,
                                                          shift = .3))),
            size = .5,
            width = .82, height = .82) +
  scale_fill_gradientn(colours = viridis::plasma(15)%>%
                         clr_desaturate(.15),
                       limits = c(0,1)) +
  guides(fill = ggfx::with_shadow(guide_colorbar(title.position = "top",
                                                 barwidth = unit(120, "pt"),
                                                 barheight = unit(5, "pt")),
    x_offset = 0, y_offset = 0, sigma = 1, color = bg_black))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        legend.position = c(.93, .05),
        legend.justification = c(1, 0),
        # legend.background = ggfx::with_shadow(element_rect(),x_offset = 0, y_offset = 0),
        legend.direction = "horizontal")

ggsave("~/Desktop/genepop.pdf",
       width = 5, height = 5, device = cairo_pdf)

