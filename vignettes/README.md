# The **vignettes** folder

This folder contains the source files of the project [documentation](https://k-hench.github.io/chapter2/).
To render the `{bookdown}` document, the working directory needs to be set to this folder.
When opening `chapter2.Rproj` with Rstudio, run:

```r
setwd('vignettes')
bookdown::render_book("index.Rmd")
```