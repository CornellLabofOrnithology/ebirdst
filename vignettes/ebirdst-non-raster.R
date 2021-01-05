## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "\\textwidth", 
                      fig.height = 5, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes local and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----load_pipd----------------------------------------------------------------
library(ebirdst)

# DOWNLOAD DATA
# Currently, example data is available on a public s3 bucket. The following 
# ebirdst_download() function copies the species results to a selected path and 
# returns the full path of the results.

# Because the non-raster data is large, there is a parameter on the
# ebirdst_download function that defaults to downloading only the raster data.
# To access the non-raster data, set tifs_only = FALSE.
sp_path <- ebirdst_download(species = "example_data", tifs_only = FALSE)
print(sp_path)

pis <- load_pis(sp_path)

