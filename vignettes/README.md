# The **vignettes** folder

This folder contains the source files of the project documentation.
To render the `{bookdown}` document, the working directory needs to be set to this folder.
When opening `chapter2.Rproj` with Rstudio, run:

```r
setwd('vignettes')
bookdown::render_book("index.Rmd")
```