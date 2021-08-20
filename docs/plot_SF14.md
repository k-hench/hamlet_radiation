---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 14



## Summary

This is the accessory documentation of Figure S14.
The Figure can be recreated by running the **R** script `plot_SF14.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF14.R 2_analysis/GxP/50000/

```

## Details of `plot_SF14.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) and on the packages [**hypoimg**](https://k-hench.github.io/hypoimg) and [**hypogen**](https://k-hench.github.io/hypogen).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF14.R 2_analysis/GxP/50000/
# ===============================================================
# This script produces Suppl. Figure 14 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/GxP/50000/')
# script_name <- "R/fig/plot_SF14.R"
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
#> ── Script: R/fig/plot_SF14.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/GxP/50000/
#> ─────────────────────────────────────────── /current/working/directory ──
```

The path containing the $G \times P$ data is received and stored inside a more descriptive variable.


```r
# config -----------------------
gxp_path <- as.character(args[1])
```

We then set the specific traits for which we want to compare the $G \times P$ association under a *linear model* (lm) versus a *linear mixed model* (lmm).


```r
# configure which gxp data to load
trait_tib  <- tibble(file = dir(gxp_path) %>% .[str_detect(.,"Bars|Peduncle|Snout")]) %>%
  mutate(prep = file) %>%
  separate(prep , into = c("trait", "model_type", "win", "step", "filetype", "zip"),
           sep = "\\.") %>%
  select(file, trait, model_type) %>%
  mutate(path = gxp_path)
```

The we load the $G \times P$ data.


```r
# load gxp data
data <- pmap_dfr(trait_tib, get_gxp_both_models)
```

After the data import, this table is already suited for plotting. 


```r
# compose final figure
p_done <- data %>%
  ggplot(aes(x = gpos, y = AVG_p_wald))+
  # add gray/white LGs background
  geom_hypo_LG()+
  # add gxp data points
  geom_point(color = plot_clr, size = .3)+
  # set axis layout
  scale_x_hypo_LG()+
  scale_fill_hypo_LG_bg()+
  # set axis titles
  labs(y = expression(G~x~P~(average~italic(p)[wald])))+
  # general plot structure separated by model type and trait
  facet_grid(trait+model_type ~ ., scales = "free_y")+
  # general plot layout
  theme_hypo()
```




Finally, we can export Figure S14.


```r
# export final figure
hypo_save(filename = "figures/SF14.png",
       plot = p_done,
       width = 11,
       height = 7,
       dpi = 600,
       type = "cairo",
       comment = plot_comment)
```

---
