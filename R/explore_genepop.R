library(GenomicOriginsScripts)
library(hypogen)
library(prismatic)
library(patchwork)

sort_element <- function(x1, x2, idx){
 sort(c(x1, x2))[[idx]] 
}

import_data <- function(file){
  read_tsv(file) %>%
    separate(Population_pair, into = c("pop1", "pop2"), sep = "-") %>%
    mutate(p = str_replace(`P-Value`, pattern = "<", replacement =  ""),
           p = str_replace(p, pattern = "Highly_sign.",
                           replacement =  "0") %>%
             as.numeric(),
           p1 = map2_chr(pop1,pop2,sort_element,idx = 1),
           p2 = map2_chr(pop1,pop2,sort_element,idx = 2))
}

plot_data <- function(dat, title){
  dat  %>%
    ggplot(aes(x = p1, y = p2, fill = p))+
    coord_equal() +
    geom_tile(aes(color = after_scale(prismatic::clr_darken(fill,shift = .3))),
              size = .5,
              width = .82, height = .82) +
    scale_fill_gradientn(colours = viridis::cividis(15)%>%
                           clr_desaturate(.15),
                         limits = c(0,1)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 barwidth = unit(120, "pt"),
                                 barheight = unit(5, "pt")))+
    labs(subtitle = title)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank(),
          legend.position = c(.93, .015),
          legend.justification = c(1, 0),
          # legend.background = ggfx::with_shadow(element_rect(),x_offset = 0, y_offset = 0),
          legend.direction = "horizontal") 
}

data_all <- import_data("2_analysis/fst_signif/workshop/genepop_summary_all.tsv")
data_ham <- import_data("2_analysis/fst_signif/workshop/genepop_summary_hamlets.tsv")
data_loc <- c("bel", "hon", "pan") %>%
  str_c("2_analysis/fst_signif/workshop/genepop_summary_", ., ".tsv") %>%
  map_dfr(import_data) %>%
  mutate(p1 = factor(p1, levels = pop_levels[-c(5,11,14)]),
         p2 = factor(p2, levels = pop_levels[-c(1,6,12)]))

# bg_black <- rgb(0, 0, 0, .75)

list(data_all, data_ham, data_loc) %>%
  map2(.y = c("all", "hamlets", "loc (only merged for ploting)"),
       plot_data) %>%
  wrap_plots(nrow = 1)

ggsave("~/Desktop/genepop.pdf",
       width = 15, height = 5, device = cairo_pdf)

# ======
slim_dat <- function(dat, col_name){  dat %>% select(p1, p2, p) %>% set_names(nm = c("p1", "p2", col_name))  }

slim_dat(data_all, "p_all") %>%
  left_join(slim_dat(data_ham, col_name = "p_ham")) %>%
  left_join(slim_dat(data_loc, col_name = "p_loc")) %>% 
  # filter(!is.na(p_loc)#,
  #        #p_ham != p_loc
  #        )
  mutate(type = ifelse(!is.na(p_loc), "loc", "non-loc")) %>%
  ggplot(aes(x = p_ham, y = p_all, color = p_loc)) +
  geom_jitter(width = .05,
              height = .05, alpha = .6)+
  scale_color_distiller(palette = "RdBu", na.value = "black")+
  facet_wrap(type ~ .)
ggsave("~/Desktop/genpop_compare.pdf", width = 8.5, height = 4, device = cairo_pdf)  



perm <- readRDS("~/Desktop/prem.rds") 
data_comb <- data_loc %>%
  mutate(run = str_c(p1, "-", p2)) %>%
  left_join(perm) %>%
  mutate(p_perm = 1 - as.double(percentile))

data_comb_summary <- data_comb %>% 
  select(run, p, p_perm) %>%
  mutate(x = ifelse(p_perm > .05, .075, .01),
         y = ifelse(p > .05, .8, .1)) %>%
  group_by(x, y) %>%
  count() %>%
  ungroup()

library(ggtext)
data_comb %>%
  ggplot(aes(x = p_perm, y = p)) +
  geom_vline(xintercept = .05) +
  geom_jitter(height = .075, alpha = .4, width = 0) +
  geom_text(data = data_comb_summary, aes(x, y, label = str_c("n: ", n))) +
  labs(y = "<i>p<sub>genepop (jittered for viz, data mostly 0 or 1)</sub></i>",
       x = "<i>p<sub>permutation (fraction of runs with fst > data)</sub></i>")+
  theme_minimal() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

ggsave("~/Desktop/fst_genepop_vs_perm.pdf", width = 6, height = 4, device = cairo_pdf)

# ================

permuation_summary <- read_rds("~/Desktop/perm_summary.rds")  %>% 
  rename(p = "p_perm") %>%
  select(p1,p2,p) %>%
  mutate(type = "permuation")


data_both_tiles <- c("bel", "hon", "pan") %>%
  str_c("2_analysis/fst_signif/workshop/genepop_summary_", ., ".tsv") %>%
  map_dfr(import_data) %>% 
  select(p1,p2,p) %>%
  mutate(type = "genepop") %>%
  bind_rows(permuation_summary) %>%
  mutate(p1 = factor(p1, levels = pop_levels),
         p2 = factor(p2, levels = pop_levels))

data_both_tiles %>%
  ggplot(aes(x = p1, y = p2)) +
  geom_tile(aes( fill = p,
                 color = after_scale(prismatic::clr_darken(fill,shift = .3))),
            size = .5,
            width = .75, height = .75) +
  geom_tile(data = data_both_tiles %>% 
              filter(p <= .05),
            fill = "transparent",
            aes(color = p),
            size = .3,
            width = .925, height = .925) +
  scale_fill_gradientn(colours = viridis::cividis(15)%>%
                         clr_desaturate(.15),
                       limits = c(0,1)) +
  scale_color_gradientn(colours = viridis::cividis(15)%>%
                          clr_desaturate(.15) %>%
                          prismatic::clr_darken(shift = .3),
                        limits =c(0, 1), guide = FALSE)+
  coord_equal()+
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = unit(120, "pt"),
                               barheight = unit(5, "pt")))+
  labs(title = "Comparison genpop vs. permutation",
       subtitle = "(upper: genepop, lower: permutation)",
       caption = "<i>p<sub>perm</sub></i> = n(<i>F<sub>ST</sub> perm</i> > <i>F<sub>ST</sub> real</i>) /1000<br>boxed tiles indicate <i>p</i> ≤ 0.05") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        legend.position = c(.93, .015),
        legend.justification = c(1, 0),
        # legend.background = ggfx::with_shadow(element_rect(),x_offset = 0, y_offset = 0),
        legend.direction = "horizontal",
        plot.caption = element_markdown()) 

ggsave("~/Desktop/fst_genepop_perm_loc.pdf", width = 5,height = 5.5, device = cairo_pdf)
