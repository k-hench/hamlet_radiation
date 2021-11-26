---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 16






## Summary

This is the accessory documentation of Figure S16.
The Figure can be recreated by running the **R** script `plot_SF16.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF16.R \
    2_analysis/admixture/ \
    metadata/phenotypes.sc
```

## Details of `plot_SF16.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF16.R \
#     2_analysis/admixture/ \
#     metadata/phenotypes.sc
# ===============================================================
# This script produces Suppl. Figure 16 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/admixture/", "metadata/phenotypes.sc")
# script_name <- "R/fig/plot_SF16.R"
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
library(paletteer)
library(patchwork)
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
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
#> ── Script: R/fig/plot_SF16.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/admixture/
#> ★ 2: metadata/phenotypes.sc
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
admx_path <- as.character(args[1])
pheno_file <- as.character(args[2])
```



```r
# load outlier window IDs (crop from admixture result file names)
gids <- dir(admx_path,
            pattern = "pop.*15.txt") %>%
  str_remove("pop.") %>%
  str_remove(".15.txt")
```



```r
# load phenotype data
pheno_data <- read_sc(pheno_file) %>%
  select(id, Bars, Peduncle, Snout) %>%
  filter(!is.na(Bars))
```



```r
# load admixture data
data <- gids %>%
  map_dfr(data_amdx, admx_path = admx_path,
          k = 2)
```



```r
# associate phenotypic trait with outlier region
pheno_facet <- tibble( trait = c("Snout","Bars",  "Peduncle"),
                       gid = c("LG04_1", "LG12_3", "LG12_4")) %>%
  mutate(facet_label = str_c(gid, " / ", trait))
```



```r
# set outlier region labels
gid_labels <-  c(LG04_1 = "LG04 (A)",
                 LG12_3 = "LG12 (B)",
                 LG12_4 = "LG12 (C)")
```



```r
# set outlier region phenotypic traits
gid_traits <-  c(LG04_1 = "Snout",
                 LG12_3 = "Bars",
                 LG12_4 = "Peduncle")
```



```r
# set path to trait images
trait_icons <- c(LG04_1 = "<img src='ressources/img/snout_c.png' width='60' />   ",
               LG12_3 = "<img src='ressources/img/bars_c.png' width='60' />    ",
               LG12_4 = "<img src='ressources/img/peduncle_c.png' width='60' />    ")
```



```r
# format phenotype data
pheno_plot_data <- data %>%
  filter(!duplicated(id)) %>%
  select(id:id_order) %>%
  left_join(pheno_data,by = c( id_nr = "id")) %>%
  arrange(spec, Bars, Peduncle, Snout, id) %>%
  mutate(ord_nr = row_number()) %>%
  pivot_longer(names_to = "trait",
               values_to = "phenotype",
               cols = Bars:Snout) %>%
  left_join(pheno_facet)
```



```r
# helper for consistent sample order across all panels
sample_order <- pheno_plot_data %>%
  filter(!duplicated(id)) %>%
  select(id, ord_nr)
```



```r
# create plot panels a-c
p_ad <- c("LG04_1", "LG12_3", "LG12_4") %>% purrr::map(adm_plot, data = data)
```



```r
# create dummy plot for the phenotype legend
p_phno <- pheno_plot_data %>%
  ggplot(aes(x = ord_nr))+
  geom_point(aes(y = trait, fill = factor(phenotype)),shape = 21)+
  scale_fill_manual("Phenotype<br><img src='ressources/img/all_traits_c.png' width='110' />",
                    values = c(`0` = "white", `1` = "black"),
                    na.value = "gray",
                    labels = c("absent", "present", "not scored"))+
  guides(fill = guide_legend(ncol = 1))+
  theme_minimal()+
  theme(legend.title = element_markdown(hjust = .5),
    legend.position = "bottom")
```



```r
# prepare table with fish annotations for the species indication
tib_drawing <- pheno_plot_data %>%
  group_by(spec) %>%
  summarise(pos = (min(ord_nr)+max(ord_nr))*.5) %>%
  ungroup()
```



```r
# create sub-plot for species indication
p_spec <- pheno_plot_data %>%
  group_by(spec) %>%
  summarise(start = min(ord_nr)-1,
            end = max(ord_nr)) %>%
  ggplot(aes(xmin = start, xmax = end,
             ymin = -Inf,
             ymax = Inf))+
  # add colored backgroud boxes
  geom_rect(aes(fill = spec), color = "black")+
  # add fish images
  (tib_drawing %>% pmap(add_spec_drawing))+
  # set axis layout
  scale_y_continuous(breaks = .5, labels = c( "Species"), limits = c(0,1))+
  scale_x_discrete(breaks = sample_order$ord_nr,
                   labels = sample_order$id,
                   expand = c(0,0)) +
  # set species color scheme
  scale_fill_manual("Species", values = clr, labels = sp_labs)+
  # set general plot layout
  theme_minimal()+
  theme(plot.title = element_text(size = 9),
        legend.position = "bottom",
        legend.text.align = 0,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())
```



```r
# create sub-plot for samplong location indication
p_loc <- pheno_plot_data %>%
  ggplot(aes(x = factor(ord_nr)))+
  # add colored boxes
  geom_raster(aes(y = 0, fill = loc))+
  # set axis layout
  scale_y_continuous(breaks = c(0),labels = c("Location"))+
  scale_x_discrete(breaks = sample_order$ord_nr,
                   labels = sample_order$id) +
  # set location color scheme
  scale_fill_manual("Location", values =  clr_loc, loc_names)+
  # set general plot layout
  theme_minimal()+
  theme(plot.title = element_text(size = 9),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())
```



```r
# compose legend from individual legend parts
p_l <- (get_legend(p_phno) %>% ggdraw()) +
  (get_legend(p_spec) %>% ggdraw()) +
  (get_legend(p_loc) %>% ggdraw()) +
  plot_layout(nrow = 1)
```



```r
# finalize figure
p_prep <- p_ad[[1]] +
  p_ad[[2]] +
  p_ad[[3]]+
  p_spec + p_loc + p_l +
  plot_layout(ncol = 1, heights = c(.4,.4,.4,.08,.02,.1)) &
  theme(legend.position = "none",
        axis.text = element_text(size = 12))
```



```r
# crop final figure (remove whitespace on left margin)
p_done <- ggdraw(p_prep, xlim = c(.023,1))
```

Finally, we can export Figure S16.


```r
# export final figure
scl <- .9
hypo_save("figures/SF16.pdf",
       plot = p_done,
       width = 16*scl,
       height = 10*scl,
       device = cairo_pdf,
       bg = "transparent",
       comment = plot_comment)
```
