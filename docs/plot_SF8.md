---
output: html_document
editor_options:
  chunk_output_type: console
---

# Supplementary Figure 8



## Summary

This is the accessory documentation of Figure S8.
The Figure can be recreated by running the **R** script `plot_SF8.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF8.R \
  2_analysis/pi/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

```

## Details of `plot_SF8.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) as well as the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**vroom**](https://vroom.r-lib.org/).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF8.R 2_analysis/pi/50k/ \
#   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script produces Suppl. Figure 8 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/pi/50k/',
#           '2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz')
# script_name <- "R/fig/plot_SF8.R"
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(vroom)
library(hypoimg)
library(hypogen)
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
#> ── Script: R/fig/plot_SF8.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/pi/50k/
#> ★ 2: 2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
#> ────────────────────────────────────────── /current/working/directory ──
```

The paths containing the $\pi$ and $\rho$ data are received and stored inside more descriptive variables.


```r
# config -----------------------
pi_path <- as.character(args[1])
rho_path <- as.character(args[2])
```

The $\pi$-path is screened for data files in the desired resolution (50 kb windows).


```r
# locate pi data files
files <- dir(pi_path, pattern = '^[a-z]{6}.50k')
```

Then, the $\pi$ data is loaded and compiled into a single table containing all populations.


```r
# load pi data
data <- str_c(pi_path, files) %>%
  purrr::map(get_pi) %>%
  bind_rows()
```

To be able to order the subplots of the final figure by the genome wide average $\pi$ of the populations, we create a second data table containing the summary of the $\pi$ data for all populations.


```r
# compute genome wide average pi for the subplot order
global_bar <- data %>%
  filter( BIN_START %% 50000 == 1) %>%
  select(N_VARIANTS, PI, spec) %>%
  group_by(spec) %>%
  summarise(genome_wide_pi = sum(N_VARIANTS*PI)/sum(N_VARIANTS)) %>%
  arrange(genome_wide_pi) %>%
  ungroup() %>%
  mutate(spec = fct_reorder(.f = spec, .x = genome_wide_pi),
         scaled_pi = genome_wide_pi/max(genome_wide_pi))
```

Then, we load the $\rho$ data.


```r
# load recombination data
rho_data <- vroom(rho_path, delim = '\t') %>%
  select(-BIN_END)
```

Next, the two data sets are merged based on the window position on the hamlet reference genome.


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

To indicate the hamlet population on the sub-plots, hamlet illustrations are loaded.


```r
# create table with fish annotations
grob_tibble2 <- global_bar$spec %>%
  purrr::map(fish_plot2) %>%
  bind_rows()
```

Then, the final figure is created.


```r
# compose final figure
p_done <- combined_data %>%
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




Finally, we can export Figure S8.


```r
# export final figure
hypo_save(filename = 'figures/SF8.pdf',
          plot = p_done,
          width = 8,
          height = 10,
          comment = plot_comment)
```

---
