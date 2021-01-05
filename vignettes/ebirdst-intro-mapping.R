## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "\\textwidth", 
                      fig.height = 4, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes locally and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----load_raster_stack--------------------------------------------------------
library(ebirdst)
library(raster)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(viridisLite)
# handle namespace conflicts
extract <- raster::extract

# Currently, example data is available on a public s3 bucket. The following 
# ebirdst_download() function copies the species results to a selected path and 
# returns the full path of the results. In this vignette, we'll use example data 
# for Yellow-bellied Sapsucker
sp_path <- ebirdst_download(species = "example_data")

# load trimmed mean abundances and label dates
abunds <- load_raster("abundance", path = sp_path)

# use parse_raster_dates() to get actual date objects for each layer
date_vector <- parse_raster_dates(abunds)
print(date_vector)

## ----aggregate----------------------------------------------------------------
abunds_low_res <- aggregate(abunds, fact = 3, fun = mean)
# native resolution
res(abunds)
# resolution after aggregation
res(abunds_low_res)

## ----project_stack------------------------------------------------------------
# define mollweide projection
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

# project single layer from stack to mollweide
week38_moll <- projectRaster(abunds[[38]], crs = mollweide, method = "ngb")
week38_moll <- trim(week38_moll)

# optionally, you can project an entire stack, but it takes much longer
# abund_moll <- projectRaster(abund, crs = mollweide, method = "ngb")

# map single layer with full annual extent
par(mar = c(0, 0, 0, 2))
plot(week38_moll, 
     col = abundance_palette(10, season = "weekly"), 
     axes = FALSE, box = FALSE, 
     maxpixels = ncell(week38_moll))

## ----map_occurrence-----------------------------------------------------------
occs <- load_raster("occurrence", path = sp_path)

# select a week in the summer
occ <- occs[[26]]

# create breaks every 0.1 from 0 to 1
occ_bins <- seq(0, 1, by = 0.1)
occ_moll <- projectRaster(occ, crs = mollweide, method = "ngb")
occ_moll <- trim(occ_moll)

par(mar = c(0, 0, 0, 2), cex = 0.9)
plot(occ_moll, 
     breaks = occ_bins, 
     col = abundance_palette(length(occ_bins) - 1, season = "weekly"),
     axes = FALSE, box = FALSE, 
     maxpixels = ncell(occ_moll),
     legend.width = 2, legend.shrink = 0.97)

## ----map_linear---------------------------------------------------------------
year_max <- max(maxValue(abunds), na.rm = TRUE)

week14_moll <- projectRaster(abunds[[14]], crs = mollweide, method = "ngb")
week14_moll <- trim(week14_moll)

# set graphical params
par(mfrow = c(1, 2), mar = c(0, 0, 0, 4))

# use raster bounding box to set the spatial extent for the plot
bb <- st_as_sfc(st_bbox(trim(week14_moll)))
plot(bb, col = "white", border = "white")
# plot the abundance
plot(week38_moll, zlim = c(0, year_max), 
     col = abundance_palette(20, season = "weekly"), 
     maxpixels = ncell(week38_moll),
     axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)

# do the same for week 14
par(mar = c(0, 0, 0, 4))
bb <- st_as_sfc(st_bbox(trim(week14_moll)))
plot(bb, col = "white", border = "white")
plot(week14_moll, zlim = c(0, year_max),
     col = abundance_palette(20, season = "weekly"), 
     maxpixels = ncell(week14_moll),
     axes = FALSE, box = FALSE, legend.shrink = 0.97, add = TRUE)

