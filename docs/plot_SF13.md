---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 13



## Summary

This is the accessory documentation of Supplementary Figure 13.
The Figure can be recreated by running the **R** script `plot_SF13.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF13.R 

```

## Details of `plot_SF13.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**ggtext**](https://wilkelab.org/ggtext/), [**hypoimg**](https://k-hench.github.io/hypoimg), [**paletteer**](https://emilhvitfeldt.github.io/paletteer/), [**patchwork**](https://patchwork.data-imaginist.com/) and [**prismatic**](https://emilhvitfeldt.github.io/prismatic/).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF13.R 2_analysis/fst_signif/random/
# ===============================================================
# This script produces Suppl. Figure 13 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/fst_signif/random/")
# script_name <- "R/fig/plot_SF13.R"
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)

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
```

```r
#> ── Script: scripts/plot_SF13.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/fst_signif/random/
#> ─────────────────────────────────────────── /current/working/directory ──
```

The directory containing the hybridization data is received and stored in a variable.


```r
# config -----------------------
rand_path <- as.character(args[1])
rand_files <- dir(path = rand_path, pattern = "_random_fst.tsv.gz")

get_random_fst <- function(file){
  nm <- file %>% str_remove(pattern = ".*\\/") %>%
    str_remove("_random_fst.tsv.gz")
  rn <- str_remove(nm, "_.*")
  sub_type <- str_remove(nm, "[a-z]{6}-[a-z]{6}_")
  vroom::vroom(file = file,
               delim = "\t",
               # col_names = c("idx", "type", "mean_fst", "weighted_fst"),
               col_types = "dcdd") %>%
    mutate(group = nm,
           run = rn,
           subset_type = sub_type)
}

data <- str_c(rand_path, rand_files) %>%
  map_dfr(.f = get_random_fst) %>% 
  mutate(sign = c(subset_non_diverged = -1, whg = 1)[subset_type])

get_percentile <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  
  # sprintf("%.7f", sum(ran < real_fst) / length(ran))
  sum(ran < real_fst) / length(ran)
}

get_n_above<- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  
  sum(ran > real_fst)
}

get_n_total <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  length(ran)
}

data_grouped <- data %>% 
  group_by(group, run, subset_type, sign) %>%
  nest() %>%
  ungroup() %>% 
  mutate(real_pop = map_dbl(data, function(data){data$weighted_fst[data$type == "real_pop"]}),
         percentile = map_chr(data, get_percentile),
         above = map_dbl(data, get_n_above),
         total = map_dbl(data, get_n_total),
         label = str_c(sprintf("%.2f", as.numeric(percentile)), "<br>(",above, "/10<sup>",log10(total),"</sup>)")) %>%
  arrange(-sign, -as.numeric(percentile), -real_pop)

run_lvls <- data_grouped %>% filter(sign == 1) %>%
  mutate(rank = row_number(),
         run = fct_reorder(run, rank)) %>% 
  .$run %>%  levels()

data_grouped_ordered <- data_grouped %>% 
  mutate(run = factor(run, levels = run_lvls))

p_done <- data %>%
  mutate(run = factor(run, levels = run_lvls)) %>%
  filter(type == "random") %>%
  ggplot(aes(group = group)) +
  geom_segment(data = data %>%
               mutate(run = factor(run, levels = run_lvls)) %>%
               filter(type == "real_pop"),
             aes(x = weighted_fst,
                 xend = weighted_fst,
                 y = 0,
                 yend = Inf * sign,
                 color = subset_type)) +
  geom_density(data = data %>%  filter(sign == 1),
                aes( x = weighted_fst, y = ..density..,
                     color = subset_type,
                     fill = after_scale(prismatic::clr_alpha(color,alpha = .3)))) +
  geom_density(data = data %>%  filter(sign == -1),
                 aes( x = weighted_fst, y = ..density.. * -1,
                      color = subset_type,
                      fill = after_scale(prismatic::clr_alpha(color,alpha = .3)))) +
  geom_richtext(data = data_grouped_ordered,
                 aes(x = .08, y = 300 * sign,
                     label = label),
                 hjust = .5, 
                 vjust = .5,
                 size = 2.5,
                 label.padding = unit(3,"pt"),
                 label.size = 0,
                 label.color = "transparent",
                 fill = "transparent") +
  scale_x_continuous(breaks = c(0, .05, .1)) +
  scale_color_manual(values = c(whg = "#4D4D4D", 
                                subset_non_diverged = "#A2A2A2"),
                     guide = FALSE) +
  labs( y = "Density (<span style='color:#4D4D4D'>whg</span> / <span style='color:#A2A2A2'>diverged regions excluded</span>)",
        x = "Weighted Average <i>F<sub>ST</sub></i>") +
  facet_wrap(run ~ ., ncol = 4, dir = "v") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

p_done <- ggdraw() + draw_plot(p_done) + 
  draw_label("Draft (10^3)", color = rgb(0,0,0,.2), size = 100, angle = 45)

scl <- 1.4
hypo_save(p_done, filename = 'figures/SF13.pdf',
          width = f_width * scl,
          height = f_width * .65 * scl,
          device = cairo_pdf,
          comment = plot_comment,
          bg = "transparent")
```
