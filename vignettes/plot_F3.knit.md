---
output: html_document
editor_options:
  chunk_output_type: console
---
# Figure 3



## Summary

This is the accessory documentation of Figure 3.
It should be possible to recreate the figure by running the **R** script `plot_F3.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_F3.R \
   2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt
```

## Details of `plot_F3.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory **R** packages [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) and on the **R** packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**vroom**](https://vroom.r-lib.org/) and [**ggforce**](https://ggforce.data-imaginist.com/).

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:



The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(ggforce)
library(hypoimg)
library(vroom)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
```

```r
#> ── Script: scripts/plot_F3.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: 2_analysis/fst/50k/
#> ★ 2: 2_analysis/summaries/fst_globals.txt
#> ─────────────────────────────────────────── /current/working/directory ──
```

The directory containing the sliding window $F_{ST}$ data and the and
the file with the genome wide average $F_{ST}$ for all the species
comparisons are received from the command line input.


```r
# config -----------------------
data_dir <- as.character(args[1])
globals_file <- as.character(args[2])
# script -----------------------
```

Then, the data folder is scanned for windowed $F_{ST}$ data with an
window size of 50 kb.


```r
# locate data files
files <- dir(path = data_dir, pattern = '.50k.windowed.weir.fst.gz')
```

Next, the genome wide average $F_{ST}$ data is loaded.


```r
# load genome wide average fst data
globals <- vroom::vroom(globals_file, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  mutate(run = str_c(loc,'-',run) %>%
           reformat_run_name()
  )
```

The package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts) 
contains the function `get_fst_fixed` to import $F_{ST}$ data and compute the
number, average length and cumulative length of regions exceeding a given $F_{ST}$
threshold.

Here, we prepare a table of a series of thresholds and all pair wise species comparisons
as a configuration table for the following import with `get_fst_fixed`.


```r
# prepare data import settings within a data table (tibble)
import_table <- list(file = str_c(data_dir,files),
                     fst_threshold = c(.5,.4,.3,.2,.1,.05,.02,.01)) %>%
  cross_df() %>%
  mutate( run =  file %>%
            str_remove('^.*/') %>%
            str_sub(.,1,11) %>%
            reformat_run_name())
```

Using the configuration table, the $F_{ST}$ data are loaded, and the
threshold-specific stats are computed.











