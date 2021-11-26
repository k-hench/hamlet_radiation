---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 9






## Summary

This is the accessory documentation of Figure S9.
The Figure can be recreated by running the **R** script `plot_SF9.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF9.R \
    2_analysis/pi/50k/ \
    2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
```

## Details of `plot_SF9.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF9.R \
#     2_analysis/pi/50k/ \
#     2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script produces Suppl. Figure 9 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/pi/50k/',
#           '2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz')
# script_name <- "R/fig/plot_SF9.R"
args <- commandArgs(trailingOnly = FALSE)
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(vroom)
library(hypoimg)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF9.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/pi/50k/
#> ★ 1: 2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
#> ────────────────────────────────────────── /current/working/directory ──
```


```r
# config -----------------------
pi_path <- as.character(args[1])
rho_path <- as.character(args[2])
```



```r
# locate pi data files
files <- dir(pi_path, pattern = '^pi.[a-z]{6}.50k')
```



```r
# load pi data
data <- str_c(pi_path, files) %>%
  purrr::map(get_pi) %>%
  bind_rows()
```



```r
# compute genome wide average pi for the subplot order
global_bar <- data %>%
  filter( BIN_START %% 50000 == 1) %>%
  select(N_SITES, PI, spec) %>%
  group_by(spec) %>%
  summarise(genome_wide_pi = sum(N_SITES*PI)/sum(N_SITES)) %>%
  arrange(genome_wide_pi) %>%
  ungroup() %>%
  mutate(spec = fct_reorder(.f = spec, .x = genome_wide_pi),
         scaled_pi = genome_wide_pi/max(genome_wide_pi))
```



```r
# load recombination data
rho_data <- vroom(rho_path, delim = '\t') %>%
  select(-BIN_END)
```



```r
# merge pi and recombination data
combined_data <- data %>%
  # filter pi data to "non-overlapping" windows
  filter(BIN_START %% 50000 == 1 ) %>%
  # reorder populations by genome wide average pi
  mutate(spec = factor(spec, levels = levels(global_bar$spec))) %>%
  # merge with recombination data
  left_join(rho_data, by = c(CHROM = 'CHROM', BIN_START = 'BIN_START'))
```



```r
# create table with fish annotations
grob_tibble2 <- global_bar$spec %>%
  purrr::map(fish_plot2) %>%
  bind_rows()
```



```r
# compose final figure
p <- combined_data %>%
  ggplot()+
  # add fish annotations
  geom_hypo_grob2(data = grob_tibble2,
                  aes(grob = grob, rel_x = .25,rel_y = .75),
                  angle = 0, height = .5,width = .5)+
  # add hex-bin desity layer
  geom_hex(bins = 30,color = rgb(0,0,0,.3),
           aes(fill=log10(..count..), x = RHO, y = PI))+
 # general plot structure (separated by run)
  facet_wrap(spec ~., ncol = 3)+
  # set axis layout and color scheme
  scale_x_continuous(name = expression(rho))+
  scale_y_continuous(name = expression(pi))+
  scico::scale_fill_scico(palette = 'berlin') +
  # customize legend
  guides(fill = guide_colorbar(direction = 'horizontal',
                               title.position = 'top',
                               barheight = unit(7,'pt'),
                               barwidth = unit(130,'pt')))+
  # general plot layout
  theme_minimal()+
  theme(legend.position = c(.84,.01),
        strip.text = element_blank())
```

Finally, we can export Figure S9.


```r
# export final figure
hypo_save(filename = 'figures/SF9.pdf',
          plot = p,
          width = 8,
          height = 10,
          comment = plot_comment)

# ===============
combined_data %>%
  filter( BIN_START %% 50000 == 1) %>%
  group_by(spec) %>%
  summarise(genom_avg_pi = sum(PI*N_SITES)/sum(N_SITES)) %>%
  write_tsv("2_analysis/summaries/pi_globals.tsv")
```
