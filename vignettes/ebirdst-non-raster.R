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

## ----map_centroids------------------------------------------------------------
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.425, 0.475))
map_centroids(sp_path, ext = lp_extent)

## ----plot_extents-------------------------------------------------------------
par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))
calc_effective_extent(sp_path, ext = lp_extent)

## ----plot_pis-----------------------------------------------------------------
# with all classes
plot_pis(pis, ext = lp_extent, by_cover_class = FALSE, n_top_pred = 15)

# aggregating fragstats for cover classes
plot_pis(pis, ext = lp_extent, by_cover_class = TRUE, n_top_pred = 15)

## ----binary_by_time-----------------------------------------------------------
plot_binary_by_time(path = sp_path, metric = "kappa", n_time_periods = 12)

## ----all_ppms-----------------------------------------------------------------
plot_all_ppms(sp_path, lp_extent)

