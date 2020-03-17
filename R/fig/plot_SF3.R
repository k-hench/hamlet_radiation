#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF2.R \
#   2_analysis/summaries/fst_globals.txt \
#   2_analysis/fst/50k/ \
#   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c( '2_analysis/summaries/fst_globals.txt',
#            '2_analysis/fst/50k/',
#            '2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(vroom)
library(furrr)
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
global_fst_file <- as.character(args[1])
fst_dir <- as.character(args[2])
rho_dir <- as.character(args[3])
# functions --------------------
# geom_hypo_grob2 <- function(mapping = NULL,
#                             data = NULL,
#                             stat = "identity",
#                             position = "identity",
#                             na.rm = FALSE,
#                             show.legend = NA,
#                             inherit.aes = FALSE,
#                             ...) {
#   layer(
#     geom = hypo_geom_grob_custom2,
#     mapping = mapping,
#     data = data,
#     stat = stat,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(na.rm = na.rm, ...)
#   )
# }
#
# hypo_geom_grob_custom2 <- ggproto(
#   "hypo_geom_grob_custom2",
#   Geom,
#   setup_data = function(self, data, params) {
#     data <- ggproto_parent(Geom, self)$setup_data(data, params)
#     data
#   },
#
#   draw_group = function(data, panel_scales, coord) {
#     vp <- grid::viewport(x=data$rel_x,
#                          y=data$rel_y,
#                          h = data$height,
#                          width = data$width,
#                          angle = data$angle)
#     g <- grid::editGrob(data$grob[[1]], vp=vp)
#     ggplot2:::ggname("geom_hypo_grob2", g)
#   },
#
#   required_aes = c("grob","rel_x","rel_y"),
#   default_aes = list(height = 1, width = 1, angle = 0)
# )
#
#
# fish_plot2 <- function(spec){
#
#   p <- ggplot()+
#     hypoimg::hypo_anno_flag(geo = loc_names[str_sub(spec,4,6)] %>% str_to_lower(),xmax = 0)+
#     hypoimg::hypo_anno_r(species = sp_names[str_sub(spec,1,3)],xmin = 0)+
#     scale_y_continuous(limits = c(-1,1))+
#     scale_x_continuous(limits = c(-.4,1))+
#     theme_void()
#
#   tibble(spec = spec, grob = list(p %>% ggplotGrob()))
# }

# load data -------------------
fst_globals <- vroom::vroom(global_fst_file, delim = '\t',
                            col_names = c('loc','run_prep','mean_fst','weighted_fst')) %>%
  separate(run_prep,into = c('pop1','pop2'),sep = '-') %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted_fst))

fst_files <- dir(fst_dir, pattern = '.50k.windowed.weir.fst.gz')

fst_data <- str_c(fst_dir,fst_files) %>%
  furrr::future_map_dfr(get_fst) %>%
  mutate(run = factor(run, levels = levels(fst_globals$run)))

# regression pi vs. rho ------------------
rho_data <- vroom::vroom(rho_dir, delim = '\t') %>%
  select(-BIN_END)

combined_data <- fst_data %>%
  filter(BIN_START %% 50000 == 1 ) %>%
  left_join(rho_data, by = c(CHROM = 'CHROM', BIN_START = 'BIN_START'))

model_data <- combined_data %>%
  group_by(run) %>%
  nest() %>%
  left_join(., fst_globals) %>%
  mutate(mod =  map(data, ~ lm(.$WEIGHTED_FST ~ .$RHO))) %>%
  bind_cols(., summarise_model(.))

grob_tibble <- model_data %>%
  ungroup() %>%
  select(pop1, pop2, loc) %>%
  set_names(c('left', 'right', 'loc')) %>%
  future_pmap_dfr(plot_fishes_location) %>%
  mutate(run = factor(run, levels = levels(fst_globals$run)))


p1 <- combined_data %>%
  ggplot()+
   geom_hypo_grob2(data = grob_tibble,
                   aes(grob = grob, rel_x = .75,rel_y = .75),
                   angle = 0, height = .5,width = .5)+
  geom_hex(bins = 30,color = rgb(0,0,0,.3),
           aes(fill=log10(..count..), x = RHO, y = WEIGHTED_FST))+
  geom_abline(data = model_data, color = rgb(1,1,1,.8),linetype = 2,
              aes(intercept = intercept, slope = slope)) +
  geom_text(data = model_data, x = 0, y = .8,parse = TRUE,hjust = 0,vjust = 1,
            aes(label = str_c('italic(R)^2:~',round(r.squared,2)))) +
  geom_text(data = model_data, x = 0, y = .95,hjust = 0,
              aes(label = run)) +
  facet_wrap(run ~., ncol = 5)+
  scale_x_continuous(name = expression(rho))+
  scale_y_continuous(name = expression(italic(F[ST])),limits = c(-.05,1))+
  scico::scale_fill_scico(palette = 'berlin') +
  guides(fill = guide_colorbar(direction = 'horizontal',
                               title.position = 'top',
                               barheight = unit(7,'pt'),
                               barwidth = unit(130,'pt')))+
  theme_minimal()+
  theme(legend.position = c(.8,.08),
        strip.text = element_blank())

p2 <- model_data %>%
  ggplot(aes(x = weighted_fst, y = slope))+
  geom_point(color = plot_clr)+
  labs(x = expression(genome~wide~weighted~mean~italic(F[ST])),
       y = expression(slope~(f(italic(F[ST]))==a~rho+b)))+
  theme_minimal()

p <- plot_grid(p1, p2,
          ncol = 1,
          rel_heights = c(1,.3),
          labels = letters[1:2] %>% project_case())

hypo_save(filename = 'figures/SF2.pdf',
          plot = p,
          width = 10,
          height = 16,
          comment = plot_comment)