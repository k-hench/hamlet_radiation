---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 6



## Summary

This is the accessory documentation of Supplementary Figure 6.
The Figure can be recreated by running the **R** script `plot_SF6.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF6.R \
  2_analysis/summaries/fst_globals.txt \
  2_analysis/fst/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

```

## Details of `plot_SF6.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) as well as the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen), [**furrr**](https://furrr.futureverse.org/), [**ggtext**](https://wilkelab.org/ggtext/) and [**vroom**](https://vroom.r-lib.org/).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF6.R \
#   2_analysis/summaries/fst_globals.txt \
#   2_analysis/fst/50k/ \
#   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script produces Suppl. Figure 6 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( '2_analysis/summaries/fst_globals.txt',
#            '2_analysis/fst/50k/',
#            '2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz')
# script_name <- "R/fig/plot_SF6.R"
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) is loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(vroom)
library(furrr)
library(ggtext)
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
#> ── Script: R/fig/plot_SF6.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/summaries/fst_globals.txt
#> ★ 2: 2_analysis/fst/50k/
#> ★ 3: 2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
#> ────────────────────────────────────────── /current/working/directory ──
```

The directories containing the $F_{ST}$ and $\rho$ data, and the file containing the genome wide $F_{ST}$ averages are
received and stored in a variable.


```r
# config -----------------------
global_fst_file <- as.character(args[1])
fst_dir <- as.character(args[2])
rho_dir <- as.character(args[3])
```

Initially, the genome wide average $F_{ST}$ values are loaded.


```r
# load genome wide average fst data
fst_globals <- vroom::vroom(global_fst_file, delim = '\t',
                            col_names = c('loc','run_prep','mean_fst','weighted_fst')) %>%
  separate(run_prep,into = c('pop1','pop2'),sep = '-') %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted_fst))
```

Then, the files containing the windowed $F_{ST}$ data are located.


```r
# locate sliding window fst data files
fst_files <- dir(fst_dir, pattern = '.50k.windowed.weir.fst.gz')
```

Next, these files are loaded.


```r
# load sliding window fst data
fst_data <- str_c(fst_dir,fst_files) %>%
  furrr::future_map_dfr(get_fst) %>%
  mutate(run = factor(run, levels = levels(fst_globals$run)))
```

After this, also the recombination data is loaded.


```r
# load recombination rate data
rho_data <- vroom::vroom(rho_dir, delim = '\t') %>%
  select(-BIN_END)
```

The $F_{ST}$ and $\rho$ data frames are merged based on the genomic position and subset to only contain non-overlapping windows.


```r
# merge fst and recombination data
combined_data <- fst_data %>%
  # filter fst data to "non-overlapping" windows
  filter(BIN_START %% 50000 == 1 ) %>%
  # merge with recombination data
  left_join(rho_data, by = c(CHROM = 'CHROM', BIN_START = 'BIN_START')) %>%
  # merge with genome wide average fst data
  left_join(.,fst_globals %>% select(run, weighted_fst)) %>%
  # add label column
  mutate(pop1 = str_sub(run,1,3),
         pop2 = str_sub(run,8,10),
         loc = str_sub(run,4,6),
         run_label = str_c("*H. ", sp_names[pop1],"* - *H. ", sp_names[pop2],"*<br>(",loc_names[loc],")" ),
         run_label = fct_reorder(run_label,weighted_fst))
```

Then, the combined data frame is nested based on the pair wise species comparison and for each comparison, $F_{ST}$ is regressed on $\rho$.
The summary statistics of these regressions are added as separate columns of the nested data frame.


```r
# nest data to run linear regression on all runs in one go
model_data <- combined_data %>%
  group_by(run) %>%
  nest() %>%
  left_join(., fst_globals) %>%
  mutate(mod =  map(data, function(data){lm(WEIGHTED_FST ~ RHO, data = data)}),
         pop1 = str_sub(run,1,3),
         pop2 = str_sub(run,8,10),
         loc = str_sub(run,4,6),
         run_label = str_c("*H. ", sp_names[pop1],"* - *H. ", sp_names[pop2],"*<br>(",loc_names[loc],")" )) %>%
  bind_cols(., summarise_model(.)) %>%
  mutate(run_label = factor(run_label, levels = levels(combined_data$run_label)))
```

Now, we can create the individual panels of the figure.

We start with panel __a__:


```r
# create subplot a (hex-bins)
p1 <- combined_data %>%
  ggplot()+
  # add hex-bin desity layer
  geom_hex(bins = 30, color = rgb(0,0,0,.3),
           aes(fill=log10(..count..),
               x = RHO, y = WEIGHTED_FST))+
  # add regression line
  geom_abline(data = model_data,
              color = rgb(1,1,1,.8),
              linetype = 2,
              aes(intercept = intercept, slope = slope)) +
  # add R^2 label
  geom_text(data = model_data, x = 0, y = .975,
            parse = TRUE, hjust = 0, vjust = 1,
            aes(label = str_c('italic(R)^2:~',round(r.squared,2)))) +
  # general plot structure (separated by run)
  facet_wrap(run_label ~., ncol = 5)+
  # set axis layout and color scheme
  scale_x_continuous(name = expression(rho))+
  scale_y_continuous(name = expression(italic(F[ST])),limits = c(-.05,1))+
  scico::scale_fill_scico(palette = 'berlin') +
  # customize legend
  guides(fill = guide_colorbar(direction = 'horizontal',
                               title.position = 'top',
                               barheight = unit(7,'pt'),
                               barwidth = unit(130,'pt')))+
  # general plot layout
  theme_minimal()+
  theme(legend.position = c(.8,.08),
        strip.text = element_markdown())
```



Then we plot the slopes of the regression within panel __b__.


```r
# create subplot b (slopes)
p2 <- model_data %>%
  ggplot()+
  geom_point(color = plot_clr,
             aes(x = weighted_fst, y = slope))+
  labs(x = expression(genome~wide~weighted~mean~italic(F[ST])),
       y = expression(slope~(f(italic(F[ST]))==a~rho+b)))+
  theme_minimal()
```



And finally the $R^2$ values in panel __c__.


```r
# create subplot c (R^2s)
p3 <- model_data %>%
  ggplot()+
  geom_point(color = plot_clr,
             aes(x = weighted_fst, y = r.squared))+
  labs(x = expression(genome~wide~weighted~mean~italic(F[ST])),
       y = expression(italic(R^2)))+
  theme_minimal()
```



The final figure is created by combining the three sub panels.


```r
# compose final figure
p_done <- plot_grid(p1,
               plot_grid(p2,p3,
                         nrow = 1,
                         labels = letters[2:3] %>%
                           project_case()),
          ncol = 1,
          rel_heights = c(1,.3),labels = project_case(c("a")))
```



Finally, we can export Figure SF6.


```r
# export final figure
hypo_save(filename = 'figures/SF6.pdf',
          plot = p_done,
          width = 10,
          height = 16,
          comment = plot_comment)
```

---
