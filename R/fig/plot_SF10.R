#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF10.R \
#     2_analysis/dxy/50k/
# ===============================================================
# This script produces Suppl. Figure 10 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/dxy/50k/")
# script_name <- "R/fig/plot_SF10.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypogen)
library(hypoimg)
library(patchwork)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
dxy_path <- as.character(args[1])

# locate dxy data files
files <- dir(dxy_path)

# load dxy data
data <- str_c(dxy_path,files) %>%
  purrr::map(get_dxy) %>%
  bind_rows() %>%
  purrr::set_names(., nm = c('scaffold', 'start', 'end', 'mid', 'sites', 'pi_pop1',
                      'pi_pop2', 'dxy', 'fst', 'GSTART', 'gpos', 'run'))

genome_wide_avg <- data %>% 
  group_by(run) %>%
  summarise(avg_dxy = mean(dxy)) %>%
  ungroup() %>%
  arrange(avg_dxy)

model_data <- data %>%
  pivot_longer(cols = starts_with("pi_pop"),
               names_to = "pi_pop", 
               values_to= "pi") %>%
  mutate(pop = str_remove(pi_pop,"pi_pop") %>% str_c("Pop. ",.))  %>%
  # filter fst data to "non-overlapping" windows
  filter(start %% 50000 == 1 ) %>%
  group_by(run, pop) %>%
  nest() %>%
  mutate(mod =  map(data, function(data){lm(pi ~ dxy, data = data)})) %>%
  bind_cols(., summarise_model(.))

dxy_subplot <- function(select_idx){
  run_select <- genome_wide_avg$run[select_idx]
  
  plt_data <- data %>%
    filter(run %in% run_select) %>%
    pivot_longer(cols = starts_with("pi_pop"),
                 names_to = "pi_pop", 
                 values_to= "pi") %>%
    mutate(pop = str_remove(pi_pop,"pi_pop") %>% 
             str_c("Pop. ",.))
  
  base_lwd <- .15
  base_line_clr <- "black"
  
  p <- plt_data %>%
    ggplot(aes(x = dxy, y = pi))+
    #coord_equal()+
    facet_grid(pop ~ run,switch = "y")+
    geom_hex(bins = 30, color = rgb(0,0,0,.3),
             aes(fill=log10(..count..)))+
    # add regression line
    geom_abline(data = model_data %>%
                  filter(run %in% run_select),
                color = rgb(1,1,1,.8),
                linetype = 2,
                aes(intercept = intercept, slope = slope)) +
    # add R^2 label
    geom_text(data = model_data%>%
                filter(run %in% run_select), x = 0, y = .022,
              parse = TRUE, hjust = 0, vjust = 1, size = 3,
              aes(label = str_c('italic(R)^2:~',round(r.squared, 3)))) +
    scale_y_continuous("\U03C0",
                       breaks = c(0,.01,.02),labels = c("0", "0.01", "0.02"))+
    scale_x_continuous(expression(italic(d[XY])),
                       breaks = c(0,.01,.02),labels = c("0", "0.01", "0.02"))+
    scico::scale_fill_scico(palette = 'berlin', limits = c(0,4.2))+
    guides(fill = guide_colorbar(direction = 'horizontal',
                                 title.position = 'top',
                                 barheight = unit(7,'pt'),
                                 barwidth = unit(130,'pt')))+
    # general plot layout
    theme_minimal()+
    theme(legend.position = "bottom",
          axis.title.y = element_text(face = "italic"),
          strip.placement = "outside", 
          strip.background.x = element_rect(fill = rgb(.95,.95,.95),
                                            colour = base_line_clr,size = base_lwd),
          panel.border = element_rect(size = base_lwd,
                                      color = base_line_clr %>% 
                                        clr_lighten(factor = .8), fill = rgb(1,1,1,0))
    )
  p
}

ps <- list(1:7, 8:14, 15:21, 22:28) %>% map(dxy_subplot)

p_done <- plot_grid(ps[[1]] + theme(legend.position = "none", axis.title.x = element_blank()), 
          ps[[2]] + theme(legend.position = "none", axis.title.x = element_blank()), 
          ps[[3]] + theme(legend.position = "none", axis.title.x = element_blank()),
          ps[[4]] + theme(legend.position = "none"),
          ps[[4]]  %>% get_legend(),
          ncol = 1, rel_heights = c(1,1,1,1,.3))

# export final figure
scl <- 1.2
hypo_save(filename = 'figures/SF10.pdf',
          plot = p_done,
          width = f_width * scl,
          height = f_width * 1.15 * scl,
          device = cairo_pdf,
          comment = plot_comment)

