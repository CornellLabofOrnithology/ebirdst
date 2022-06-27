
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ebirdst: Access and Analyze eBird Status and Trends Data

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/CornellLabofOrnithology/ebirdst/workflows/R-CMD-check/badge.svg)](https://github.com/CornellLabofOrnithology/ebirdst/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/ebirdst)](https://cran.r-project.org/package=ebirdst)
<!-- badges: end -->

## Installation

Install `ebirdst` from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("ebird/ebirdst", build = FALSE)
```

This version of `ebirdst` is designed to work with the eBird Status Data
Products estimated for the year 2020, with visualizations being released
on the web in November 2021, and data access being made available in
June 2022 **Users are strongly discouraged from comparing Status and
Trends results between years due to methodological differences between
versions.** If you have accessed and used previous versions and/or may
need access to previous versions for reasons related to reproducibility,
please contact <ebird@cornell.edu> and your request will be considered.

## Data access

Data access is granted through an Access Request Form at:
<https://ebird.org/st/request>. Access with this form generates a key to
be used with this R package and is provided immediately (as long as
commercial use is not requested). Our terms of use have been designed to
be quite permissive in many cases, particularly academic and research
use. When requesting data access, please be sure to carefully read the
terms of use and ensure that your intended use is not restricted.

After completing the Access Request Form, you will be provided a Status
and Trends access key, which you will need when downloading data. To
store the key so the package can access it when downloading data, use
the function `set_ebirdst_access_key("XXXXX")`, where `"XXXXX"` is the
access key provided to you. **Restart R after setting the access key.**

## Citation

If you use the the eBird Status & Trends data please cite it with:

<blockquote>
Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S.
Ligocki, W. Hochachka, L. Jaromczyk, C. Wood, I. Davies, M. Iliff, L.
Seitz. 2021. eBird Status and Trends, Data Version: 2020; Released:
2021. Cornell Lab of Ornithology, Ithaca, New York.
<a href="https://doi.org/10.2173/ebirdst.2020" class="uri">https://doi.org/10.2173/ebirdst.2020</a>
</blockquote>

## Vignettes

For full package documentation, including a series of vignettes covering
the full spectrum from introductory to advanced usage, please see the
package [website](https://cornelllabofornithology.github.io/ebirdst).
The available vignettes are:

-   [Introduction to eBird Status & Trends
    Data](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst.html):
    covers data access, available data products, and structure and
    format of data files.
-   [Working with Raster
    Data](https://cornelllabofornithology.github.io/ebirdst/articles/rasters.html):
    loading and analysing the raster data products.

## Quick Start

This quick start guide shows how to download data and plot abundance
values similar to how they are plotted for the [eBird Status and Trends
weekly abundance
animations](https://ebird.org/science/status-and-trends/yebsap/abundance-map-weekly).
In this guide, and throughout all package documentation, a simplified
example dataset is used consisting of Yellow-bellied Sapsucker in
Michigan. For a full list of the species available for download, look at
the data frame `ebirst_runs`, which is included in this package.

**Important note: after downloading the results, do not change the file
structure.** All functionality in this package relies on the structure
inherent in the delivered results. Changing the folder and file
structure will cause errors with this package.

``` r
library(ebirdst)
library(raster)
library(sf)
library(fields)
library(rnaturalearth)

# download example data, yellow-bellied sapsucker in michigan
path <- ebirdst_download(species = "example_data")

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster(path = path, resolution = "lr")

# load species specific mapping parameters
pars <- load_fac_map_parameters(path)
# custom coordinate reference system
crs <- pars$custom_projection
# legend breaks
breaks <- pars$weekly_bins
# legend labels for top, middle, and bottom
labels <- pars$weekly_labels

# get a date vector specifying which week each raster layer corresponds to
weeks <- parse_raster_dates(abd)
print(weeks)
#>  [1] "2020-01-04" "2020-01-11" "2020-01-18" "2020-01-25" "2020-02-01" "2020-02-08"
#>  [7] "2020-02-15" "2020-02-22" "2020-03-01" "2020-03-08" "2020-03-15" "2020-03-22"
#> [13] "2020-03-29" "2020-04-05" "2020-04-12" "2020-04-19" "2020-04-26" "2020-05-03"
#> [19] "2020-05-10" "2020-05-17" "2020-05-24" "2020-05-31" "2020-06-07" "2020-06-14"
#> [25] "2020-06-21" "2020-06-28" "2020-07-06" "2020-07-13" "2020-07-20" "2020-07-27"
#> [31] "2020-08-03" "2020-08-10" "2020-08-17" "2020-08-24" "2020-08-31" "2020-09-07"
#> [37] "2020-09-14" "2020-09-21" "2020-09-28" "2020-10-05" "2020-10-12" "2020-10-19"
#> [43] "2020-10-26" "2020-11-02" "2020-11-09" "2020-11-16" "2020-11-23" "2020-11-30"
#> [49] "2020-12-07" "2020-12-14" "2020-12-21" "2020-12-28"

# select a week in the middle of the year
abd <- abd[[26]]

# project to species specific coordinates
# the nearest neighbor method preserves cell values across projections
abd_prj <- projectRaster(abd, crs = crs, method = "ngb")

# get reference data from the rnaturalearth package
# the example data currently shows only the US state of Michigan
wh_states <- ne_states(country = c("United States of America", "Canada"),
                       returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# start plotting
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

# use raster bounding box to set the spatial extent for the plot
bb <- st_as_sfc(st_bbox(trim(abd_prj)))
plot(bb, col = "white", border = "white")
# add background reference data
plot(wh_states, col = "#cfcfcf", border = NA, add = TRUE)

# plot zeroes as light gray
plot(abd_prj, col = "#e6e6e6", maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# define color palette
pal <- abundance_palette(length(breaks) - 1, "weekly")
# plot abundance
plot(abd_prj, col = pal, breaks = breaks, maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# state boundaries
plot(wh_states, add = TRUE, col = NA, border = "white", lwd = 1.5)

# legend
label_breaks <- seq(0, 1, length.out = length(breaks))
image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
           smallplot = c(0.90, 0.93, 0.15, 0.85),
           legend.only = TRUE,
           axis.args = list(at = c(0, 0.5, 1), 
                            labels = round(labels, 2),
                            cex.axis = 0.9, lwd.ticks = 0))
```

<img src="man/figures/README-quick_start-1.png" width="100%" />
