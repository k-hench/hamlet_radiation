library(GenomicOriginsScripts)
library(prismatic)
library(ggtext)

legend_bivariate <- function(x = c(1, 15),
                             y = c(1, 15), 
                             x_title = "x",
                             y_title = "y",
                             dir_x = 1,
                             dir_y = 1,
                             pal = clr,
                             fun, n_bins = 15){
  data <- expand.grid(x_d = 1:n_bins,
                      y_d = 1:n_bins)
  
  x_r <- range(x, na.rm = TRUE)
  y_r <- range(y, na.rm = TRUE)
  x_lab <- x_r %>% scales::breaks_pretty()(3)
  y_lab <- y_r %>% scales::breaks_pretty()(3)
  x_br <- scales::rescale(x_lab, from = x_r, to = c(1,n_bins))
  y_br <- scales::rescale(y_lab, from = y_r, to = c(1,n_bins))
  x_f <- dplyr::between(x_br, left = 1, right = n_bins)
  y_f <- dplyr::between(y_br, left = 1, right = n_bins)
  
  data %>% 
    ggplot(aes(x = x_d, y = y_d )) +
    geom_raster(aes(fill = x_d, alpha = y_d* dir_y)) +
    scale_fill_gradientn(colours = pal) +
    scale_x_continuous(name = x_title, breaks = x_br, labels = x_lab,
                       position = "top") +
    scale_y_continuous(name = y_title, breaks = y_br, labels = y_lab) +
    coord_equal(expand = 0) +
    # theme_void() +
    theme_minimal(base_size = plot_text_size) +
    theme(axis.title.x = element_markdown(margin = margin(5,5,5,5)),
          axis.title.y = element_markdown(margin = margin(5,5,5,5)),
          axis.line = element_line(size = .2),
          axis.ticks = element_line(size = .2),
          panel.grid = element_blank(),
          legend.position = "none")
}

two_chr_to_sorted_pair <- function(P2, P3, ...){
  n1 <- as.numeric(factor(P2, levels = sorter$group))
  n2 <- as.numeric(factor(P3, levels = sorter$group))
  char_sorted <- if(n1 < n2){c(P2, P3)} else { c(P3, P2)}
  str_c(char_sorted[[1]], "-", char_sorted[[2]])
}

sorter <- read_tsv("ressources/species_order_tree1.txt",
                   col_names = "group")

p_cap <- 8

data <- read_tsv("2_analysis/dstats/hyp_ld05_dtrios_BBAA.txt") %>%
    mutate(pair = pmap_chr(.,two_chr_to_sorted_pair)) %>% 
    group_by(pair) %>% 
    mutate(p_is_min = `p-value` == min(`p-value`)) %>%
    filter(p_is_min) %>% 
    mutate(d_is_max = Dstatistic == max(Dstatistic)) %>% 
    filter(d_is_max) %>% 
    ungroup() %>% 
    separate(pair, into = c("p_left", "p_right"), sep = "-", remove = FALSE) %>% 
  mutate(`p-value` = if_else(`p-value` < 10^-p_cap,10^-p_cap,`p-value` ))

data_prep <- data %>%
  dplyr::select(p_left, p_right, Dstatistic, `p-value`) %>% 
  bind_rows(data %>%
              dplyr::select(p_left = p_right, p_right = p_left,
                            Dstatistic, `p-value`)) 

data_full <- cross_df(list(p_left = sorter$group , 
                           p_right = sorter$group)) %>%
  left_join(data_prep) %>% 
  mutate(p_left = factor(p_left, levels = sorter$group),
         p_right = factor(p_right, levels = rev(sorter$group)),
         Dstatistic = if_else(is.na(Dstatistic), 0, Dstatistic),
         `p-value` = if_else(is.na(`p-value`), 1, `p-value`))

data_signif <- read_tsv("2_analysis/dstats/BBAA_sign_ld05.csv") %>%
  mutate(pair = pmap_chr(.,two_chr_to_sorted_pair)) %>% 
  separate(pair, into = c("p_left", "p_right"), sep = "-", remove = FALSE) %>% 
  mutate(p_left = factor(p_left, levels = sorter$group),
         p_right = factor(p_right, levels = rev(sorter$group)))

clr <- scales::colour_ramp(colors = c("#02CAEE", "#7D6181", "#EC041F"))((1:7)/7) %>% clr_saturate(.1)
clr <- scales::colour_ramp(colors = RColorBrewer::brewer.pal(5,"PuOr")[c(1,2,4,5)])((1:7)/7) %>% clr_saturate(.1)
# clr <- viridis::cividis(6)
# clr <- scico::scico(5,palette = "tokyo")

d_lim <- c(0,.01)
p_lim <- c(1, p_cap)
p2 <- legend_bivariate(x = d_lim,
                       y = -p_lim,
                       fun = clr_lighten,
                       pal = clr,
                       dir_y = -1,
                       y_title = "log<sub>10</sub>*(p)*",
                       x_title = "D",
                       n_bins = 23)

(p1 <- data_full %>%
  mutate(check1 = as.numeric(p_left),
         check2 = as.numeric(p_right)) %>% 
  filter(check2 < 17-check1) %>% 
  ggplot()+
  geom_tile(aes(x = p_left, 
                y = p_right, 
                fill = Dstatistic,
                alpha = -log10(`p-value`),
                color = after_scale(clr_darken(fill)))) +#log10(`p-value`))) +
    # geom_text(data = data_signif %>% group_by(p_left,p_right) %>% count(),
    #             aes(x = p_left,
    #                 y = p_right,
    #                 label = n),
    #           size = plot_text_size / .pt,
    #           color = "white") +
    geom_jitter(data = data_signif,
              width = .2, height = .2,
              aes(x = p_left,
                  y = p_right,
                  size = -holm,
                  fill = Dstatistic,
                  color = after_scale(clr_darken(fill))),
              shape = 21, alpha = .5) +
  annotation_custom(grob = ggplotGrob(p2),xmin = 8,ymin = 8) +
  scale_fill_gradientn(colours = clr,
                       limits = d_lim,
                       guide = "none"
                       ) +
  scale_alpha(limits = p_lim,
              na.value = 0, guide = "none"
              ) +
  scale_size(range = c(.1, 2.5),
             guide =  "none") +
  coord_equal() +
  theme_minimal(base_size = plot_text_size) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5),
        axis.title = element_blank()) )

ggsave(filename = "figures/SFx5.pdf",
       width = f_width_half,
       height = f_width_half,
       device = cairo_pdf, bg = "transparent")

