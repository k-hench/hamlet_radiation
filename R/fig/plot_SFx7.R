#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFx7.R figures/data/whg_tree/
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c('figures/data/whg_tree/')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(ggtree)
library(ape)
library(phytools)
library(hypoimg)
library(geomfactory)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)

# config -----------------------
tree_path <- as.character(args[1])
geomfactory::factory_geom_point('support')

# functions -----------------

tree_select <- get_tree(loc = 'all',mode = 'whg_no_og', tree_dir = tree_path)$tree[[1]]

spec_list <- list('uni','pue','may','nig','gem')
clr_plt <- c(clr2,gem = "#414D6BFF",ungrouped = "darkgray")
clr_plt['may'] <- 'darkblue'

tree_rot <- ggtree(tree_select, layout = 'circular') %>%
  ggtree::rotate(463) %>%
  ggtree::rotate(551) %>%
  ggtree::rotate(653) %>%
  ggtree::rotate(397) %>%
  ggtree::rotate(356) %>%
  ggtree::rotate(415) %>%
  ggtree::rotate(419) %>%
  ggtree::rotate(663) %>%
  .$data

tree_df <- tree_rot %>%
  mutate(id = str_remove(label,"_[A,B]"),
         spec = str_sub(id,start = -6,end = -4) %>% str_remove(.,"[0-9.]{1,3}$")%>% str_remove(.," "),
         spec = ifelse(spec == "",'ungrouped',spec),
         geo = str_sub(id,start = -3,end = -1) %>% str_remove(.,"[0-9]{1,3}$")%>% str_remove(.," "),
         geo = ifelse(geo == "",'ungrouped',geo),
         grouped = ifelse(spec == 'ungrouped', 'ungrouped', 'species'))

lyout <- 'circular'

geo_chunks <- tree_df %>%
  filter(isTip) %>%
  select(y,spec,geo) %>%
  arrange(y) %>%
  mutate(check = 1-(geo == lag(geo,default = '')),
         chunk = cumsum(check)) %>%
  group_by(chunk) %>%
  summarise(geo = geo[1],
            ymin = min(y),
            ymax = max(y))

c_vals <- c(ungrouped = rgb(.6, .6, .6),
            R = rgb(.6, .6, .6), clr2)
c_labs <- c(R = "",ungrouped = "",sp_labs)
c_breaks <- c(NA, NA, names(clr2))

p <- ggtree(tree_df, layout = lyout,
         aes(color = spec), color = NA)+
  geom_rect(inherit.aes = FALSE,
            data = geo_chunks,
            aes(xmin = .0035, xmax = .35,
                ymin = ymin, ymax = ymax,
                fill = geo), size = 2)+
  geom_rect(inherit.aes = FALSE,
            data = tibble(xmax = .335),
            aes(xmin = -Inf, xmax = xmax,
                ymin = -Inf, ymax = Inf),
            fill = rgb(1, 1 ,1 ,.9), size = 2)+
  geom_tree(data = tree_df,
            layout=lyout) +
  geom_tippoint(size = .4)+
  geom_nodepoint_support(data = tree_df %>% filter(!isTip),
                         aes(support_f = as.numeric(label)),
                         size = 1,
                         shape = 21)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,.5))+
  scale_color_manual(values = c_vals,
                     breaks = c_breaks,
                     labels = c_labs)+
  scale_fill_manual(values = clr_loc,
                    breaks = names(clr_loc),
                    labels = loc_names)+
  scale_support_f_continuous(low = 'lightgray',
                             high = 'black',
                             guide = guide_colourbar_support(title = 'Node support',
                                                             direction = 'horizontal',
                                                             order = 3,
                                                             title.position = 'top',
                                                             barheight = unit(5, 'pt'),
                                                             barwidth = unit(150, 'pt')))+
  guides(fill = guide_legend(title = 'Location', ncol = 2, order = 2),
         color = guide_legend(title = 'Species',
                              ncol = 2,
                              order = 1,
                              label.hjust = 0))+
  theme_void()

hypo_save(filename = 'figures/SX7.pdf',
          plot = p,
          width = 9,
          height = 6,
          comment = plot_comment)
