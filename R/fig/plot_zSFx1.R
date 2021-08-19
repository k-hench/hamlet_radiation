library(tidyverse)
library(GenomicOriginsScripts)
library(ggtext)
svardal <- readxl::read_xlsx("ressources/Svardal_etal_2021.xlsx")
leffler <- read_tsv("ressources/Leffler_etal_S1.txt")

# spec_rates <- read_csv("~/Downloads/pre_icrs/FToL_Rates_Speciation_top100.csv")
# 
# spec_rates %>% 
#   ggplot(aes(x = `lambda.tc`)) +
#   geom_histogram()

# ===========
leffler %>% 
  select( `Latin name`:`Common name`,
          `Mean diversity/bp across loci (%); (populations separated by semi-colon)`,
          `ThetaW or Pi`) %>% 
  # filter(`Common name` %in% c("Eurasian lynx","Human",  "Three-spined stickleback",
  #                             "Collared flycatcher"))
  .$`Common name` %>%  unique() %>% sort()


pi_globals <- read_tsv("2_analysis/summaries/pi_globals.tsv") %>% 
  dplyr::select(Species = "spec" , Pi = "genome_wide_pi") %>% 
  mutate(`Radiation/group` = sp_names[str_sub(Species,1,3)])

hypo_lvls <- pi_globals$`Radiation/group` %>% unique()

lvls <-  c("Victoria", "Malawi", "Mweru", "Other", "Tanganyika",
           "Natron", "Nicaragua", "ParanaUruguay", "Hypoplectrus")

hypo_groups <- tibble(`Radiation/group` = hypo_lvls,
       group = "Hypoplectrus",
       group_f = factor(group, lvls),
       name = str_c("*H. ",`Radiation/group`,"*"),
       x = 1 - seq_along(hypo_lvls))

svardal_gr %>% tail(18)
svardal_gr <- tibble(`Radiation/group` = svardal$`Radiation/group` %>%  unique(),
       group = c("Tanganyika", "Victoria", "Tanganyika", "Other", "Victoria",
                 "Tanganyika","Tanganyika", "Victoria", "Malawi", "Malawi",
                 "Tanganyika", "Natron", "Nicaragua", "Tanganyika", "Natron",
                 "Tanganyika", "Mweru", "ParanaUruguay", "Tanganyika", "Nicaragua",
                 "Tanganyika", "Malawi", "Tanganyika", "Natron", "Tanganyika",
                 "Tanganyika", "Natron"),
       group_f = factor(group, lvls),
       name = str_remove(`Radiation/group`, pattern = ".*-") %>% 
         str_remove(.,"^ ") %>% 
         ifelse(. == "A. calliptera", "*A. calliptera*", .),
       x = c(17,3,14,8,2, 9,19,1,4,5,
             13,23,26,12,22,15,7,27,18,25,
             16,6,10,24,20, 11,21)
       ) %>% 
  bind_rows(hypo_groups)

svardal_combined <- svardal %>% 
  select(Species, `Radiation/group`, Pi) %>%
  bind_rows(pi_globals) %>%
  left_join(svardal_gr)

svardal_summary <- svardal_combined %>% 
  group_by(`Radiation/group`) %>% 
  summarise(pi_min = min(Pi),
            pi_max = max(Pi),
            pi_mean = mean(Pi),
            n = length(Pi)) %>% 
  ungroup() %>% 
  left_join(svardal_gr)

group_summary <- svardal_combined %>% 
  group_by(`Radiation/group`) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  group_by(group) %>% 
  summarise(x_min = min(x),
            x_max = max(x),
            x_mean = mean(x)) %>% 
  ungroup() %>% 
  mutate(angle = ifelse(group %in% c(#"Hypoplectrus",
                                     "Mweru", "Other",
                                     "Nicaragua", "ParanaUruguay"), 90, 0))

svardal_combined %>% 
  ggplot(aes(x = x,
             y = Pi)) +
  geom_rect(data = group_summary,
            inherit.aes = FALSE,
            aes(xmin = x_min - .4, xmax = x_max +.4,
                ymin = -Inf, ymax = Inf),
            fill = rgb(.7,.7,.7,.3)) +
  geom_point(data = svardal_summary %>% filter(n > 1), 
             aes(y = pi_mean),
             size = 2, color ="red") +
  geom_segment(data = svardal_summary %>% filter(n > 1), 
               aes(xend = x,
                   y = pi_min, yend = pi_max)) +
  geom_point(size = .8) +
  geom_richtext(data = svardal_summary,
            aes(x = x, y = pi_max + .0001, label = name),
            size = 3, angle = 90, hjust = 0, color = rgb(0,0,0,.5),
            fill = "transparent", label.size = 0, label.color = "transparent") +
  geom_text(data = group_summary %>%  filter(angle == 0),
            aes(x = x_mean, label = group, y = Inf),
            vjust = 1.8, hjust = .5, fontface = "bold", color = rgb(0,0,0,.5)) +
  geom_text(data = group_summary %>%  filter(angle == 90),
            aes(x = x_mean, label = group, y = Inf),
            angle = 90,
            vjust = .5, hjust = 1.1, fontface = "bold", color = rgb(0,0,0,.5)) +
  scale_x_continuous(breaks = svardal_summary$x, labels = svardal_summary$name) +
  # scale_fill_brewer(palette = "Set1", guide = FALSE) +
  # scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "gray"), guide = FALSE) +
  scale_y_continuous(name = "\U03C0") +
  coord_cartesian(ylim = c(0,.0078),expand = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave("figures/SFx1.pdf",width = 12,height = 4, device = cairo_pdf, bg = "transparent")
