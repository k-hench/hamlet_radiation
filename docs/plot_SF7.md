---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 7






## Summary

This is the accessory documentation of Figure S7.
The Figure can be recreated by running the **R** script `plot_SF7.R`:

```sh
cd $BASE_DIR

run from terminal:
Rscript --vanilla R/fig/plot_SF7.R \
    2_analysis/dxy/50k/
```

## Details of `plot_SF7.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF7.R \
#     2_analysis/dxy/50k/
# ===============================================================
# This script produces Suppl. Figure 7 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/dxy/50k/')
# script_name <- "R/fig/plot_SF7.R"
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
library(hypoimg)
library(hypogen)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
```

```r
#> ── Script: R/fig/plot_SF7.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/dxy/50k/
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
dxy_path <- as.character(args[1])
```



```r
# locate dxy data files
files <- dir(dxy_path)
```



```r
# load dxy data
data <- str_c(dxy_path,files) %>%
  purrr::map(get_dxy) %>%
  bind_rows() %>%
  purrr::set_names(., nm = c('scaffold', 'start', 'end', 'mid', 'sites', 'pi_pop1',
                      'pi_pop2', 'dxy', 'fst', 'GSTART', 'gpos', 'run'))
```



```r
# create table for the indication of genome wide average dxy in the plot background
# (rescale covered dxy range to the extent of the genome)
global_bar <- data %>%
  # filter to non-overlaping windows only
  filter( start %% 50000 == 1) %>%
  select(sites, dxy, run) %>%
  group_by(run) %>%
  summarise(genome_wide_dxy = sum(sites*dxy)/sum(sites)) %>%
  arrange(genome_wide_dxy) %>%
  ungroup() %>%
  mutate(run = fct_reorder(.f = run, .x = genome_wide_dxy),
         scaled_dxy = genome_wide_dxy/max(genome_wide_dxy))
```



```r
# prepare plot annotaton images
grob_tibble <-  global_bar %>%
  mutate(loc = str_sub(run,4,6),
         right = str_sub(run,1,3),
         left = str_sub(run,8,10)) %>%
  select(1,4:6) %>%
  pmap(.,plot_pair_run) %>%
  bind_rows()
```



```r
# prepare plotting elements --------
# pre-define secondary x-axis breaks
sc_ax <- scales::cbreaks(c(0,max(global_bar$genome_wide_dxy)),
                         scales::pretty_breaks(4))
```



```r
# pre-define secondary x-axis labels
labels <- str_c(c("", sc_ax$breaks[2:5]*1000),
                c("0", rep("\u00B710^-3",4)))
```



```r
# sort pair-wise population comparisons by average genome wide dxy
data <- data %>%
  mutate(run = factor(run, levels = levels(global_bar$run)))
```



```r
# compose final figure
p_done <- ggplot()+
  # general plot structure separated by run
  facet_wrap( .~run, as.table = TRUE, ncol = 1, dir = 'v')+
  # add genome wide average dxy in the background
  geom_rect(data = global_bar %>% mutate(xmax = scaled_dxy * hypo_karyotype$GEND[24]),
            aes(xmin = 0, xmax = xmax, ymin = -Inf, ymax = Inf), color = rgb(1,1,1,0),fill = clr_below)+
  # add LG borders
  geom_vline(data = hypogen::hypo_karyotype, aes(xintercept = GEND), color = hypo_clr_lg)+
  # add dxy data points
  geom_point(data = data, aes(x = gpos, y = dxy),
             size=.2,color = plot_clr) +
  # add fish images
  geom_hypo_grob2(data = grob_tibble,
                  aes(grob = grob, rel_x = .945, rel_y = .5),
                  angle = 0, height = .9, width = .13)+
  # axis layout
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[24],
                                       breaks = (sc_ax$breaks/max(global_bar$genome_wide_dxy))[1:5],
                                       labels = labels,
                                       name = expression(Genomic~position/~Genome~wide~italic(d[XY]))))+
  scale_y_continuous(name = expression(italic(d[XY])), breaks = c(0,.01, .02))+
  # set plot extent
  coord_cartesian(xlim = c(0, hypo_karyotype$GEND[24]*1.135))+
  # general plot layout
  theme_hypo()+
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(),
        axis.text.x.bottom = element_markdown(colour = 'darkgray'))
```

Finally, we can export Figure S7.


```r
# export final figure
hypo_save(filename = 'figures/SF7.png',
          plot = p_done,
          width = 8,
          height = 12,
          dpi = 600,
          type = "cairo",
          comment = plot_comment)

system("convert figures/SF7.png figures/SF7.pdf")
system("rm figures/SF7.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SF7.pdf")
system(create_metadata)
```
