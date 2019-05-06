<!-- README.md is generated from README.Rmd. Please edit that file -->
ebirdst: Access and Analyze eBird Status and Trends Data
========================================================

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![AppVeyor Build
status](https://ci.appveyor.com/api/projects/status/v7cyxwquwrxxa5l6/branch/master?svg=true)](https://ci.appveyor.com/project/mstrimas/ebirdst/branch/master)
[![Travis build
status](https://travis-ci.org/CornellLabofOrnithology/ebirdst.svg?branch=master)](https://travis-ci.org/CornellLabofOrnithology/ebirdst)
[![Coverage
status](https://codecov.io/gh/CornellLabofOrnithology/ebirdst/branch/master/graph/badge.svg)](https://codecov.io/github/CornellLabofOrnithology/ebirdst?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/ebirdst)](https://cran.r-project.org/package=ebirdst)

Installation
------------

Install `ebirdst` from CRAN with:

    install.packages("ebirdst")

Alternatively, you can install the development version from GitHub with:

    # install.packages("remotes")
    remotes::install_github("CornellLabofOrnithology/ebirdst")

Vignettes
---------

For a full introduction and advanced usage, please see the package
[website](https://cornelllabofornithology.github.io/ebirdst). An
[introductory
vignette](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-introduction.html)
details the data access and structure of the results. An [intro mapping
vignette](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-intro-mapping.html)
expands upon the quick start readme and shows the basic mapping moves.
The [advanced mapping
vignette](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-advanced-mapping.html)
shows how to reproduce the seasonal maps and statistics on the [eBird
Status and Trends website](https://ebird.org/science/status-and-trends).
Finally, the [non-raster data
vignette](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-pipd.html)
details how to access additional information from the model results
about predictor importance and directionality, as well as predictive
performance metrics.

Quick Start
-----------

This quick start guide shows how to download data and plot abundance
values similar to how they are plotted for the [eBird Status and Trends
Abundance
animations](https://ebird.org/science/status-and-trends/woothr/abundance-map-weekly).
In this guide, a simplified example dataset is used consisting of
Yellow-bellied Sapsucker in Michigan. For a full list of the species
available for download, look at the data frame `ebirst_runs`, which is
included in this package.

**Important note: after downloading the results, do not change the file
structure.** All functionality in this package relies on the structure
inherent in the delivered results. Changing the folder and file
structure will cause errors with this package. If you use this package
to analyze the results, you do not ever need to interact with the files
directly, outside of R.

    library(ebirdst)
    library(viridis)
    library(raster)
    library(sf)
    library(rnaturalearth)

    # download data
    # download a simplified example dataset from aws s3
    # example data are for yellow-bellied sapsucker in michigan
    # fby default ile will be stored in a persistent data directory:
    # rappdirs::user_data_dir("ebirdst"))
    sp_path <- ebirdst_download(species = "example_data")

    # load estimated relative abundance and label with dates
    # this raster stack has 52 layers, one for each week of the year
    abunds <- load_raster("abundance_umean", path = sp_path)

    # use parse_raster_dates() to get actual date objects for each layer
    date_vector <- parse_raster_dates(abunds)
    print(date_vector)
    #>  [1] "2016-01-04" "2016-01-11" "2016-01-18" "2016-01-25" "2016-02-01"
    #>  [6] "2016-02-08" "2016-02-15" "2016-02-22" "2016-03-01" "2016-03-08"
    #> [11] "2016-03-15" "2016-03-22" "2016-03-29" "2016-04-05" "2016-04-12"
    #> [16] "2016-04-19" "2016-04-26" "2016-05-03" "2016-05-10" "2016-05-17"
    #> [21] "2016-05-24" "2016-05-31" "2016-06-07" "2016-06-14" "2016-06-21"
    #> [26] "2016-06-28" "2016-07-06" "2016-07-13" "2016-07-20" "2016-07-27"
    #> [31] "2016-08-03" "2016-08-10" "2016-08-17" "2016-08-24" "2016-08-31"
    #> [36] "2016-09-07" "2016-09-14" "2016-09-21" "2016-09-28" "2016-10-05"
    #> [41] "2016-10-12" "2016-10-19" "2016-10-26" "2016-11-02" "2016-11-09"
    #> [46] "2016-11-16" "2016-11-23" "2016-11-30" "2016-12-07" "2016-12-14"
    #> [51] "2016-12-21" "2016-12-28"

    # select a week in the summer
    abund <- abunds[[26]]
    rm(abunds)

    # project to mollweide for mapping
    # the nearest neighbor method preserves cell values across projections
    mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
    abund_moll <- projectRaster(abund, crs = mollweide, method = "ngb")

    # get reference data from the rnaturalearth package
    # the example data currently shows only the US state of Michigan
    wh_states <- ne_states(country = c("United States of America", "Canada"),
                        returnclass = "sf") %>% 
      st_transform(crs = mollweide) %>% 
      st_geometry()
      
    # calculate ideal color bins for abundance values for this week
    week_bins <- calc_bins(abund_moll)

    # start plotting
    par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))

    # use raster bounding box to set the spatial extent for the plot
    bb <- st_as_sfc(st_bbox(trim(abund_moll)))
    plot(bb, col = "white", border = "white")

    # add background reference data
    plot(wh_states, col = "#eeeeee", border = NA, add = TRUE)

    # plot zeroes as gray
    plot(abund_moll == 0, col = "#dddddd", 
         maxpixels = ncell(abund_moll),
         axes = FALSE, legend = FALSE, add = TRUE)

    # define color bins
    qcol <- abundance_palette(length(week_bins$bins) - 1, "weekly")

    # plot abundances
    plot(abund_moll, col = qcol, breaks = week_bins$bins,
         maxpixels = ncell(abund_moll),
         axes = FALSE, legend = FALSE, add = TRUE)

    # for legend, create a smaller set of bin labels
    bin_labels <- format(round(week_bins$bins, 2), nsmall = 2)
    bin_labels[!(bin_labels %in% c(bin_labels[1],
                                   bin_labels[round((length(bin_labels) / 2)) + 1],
                                   bin_labels[length(bin_labels)]))] <- ""
    bin_labels <- c("0", bin_labels)

    # create colors that include gray for 0
    lcol <- c("#dddddd", qcol)

    # set legend such that color ramp appears linearly
    ltq <- seq(from = week_bins$bins[1], to = week_bins$bins[length(week_bins$bins)],
               length.out = length(week_bins$bins))
    ltq <- c(0, ltq)

    # plot legend
    plot(abund_moll ^ week_bins$power, legend.only = TRUE,
         col = lcol,
         breaks = ltq ^ week_bins$power, 
         lab.breaks = bin_labels, legend.shrink = 0.97,
         legend.width = 2, axis.args = list(cex.axis = 0.9, lwd.ticks = 0))

    # add state boundaries on top
    plot(st_geometry(wh_states), add = TRUE, col = NA, border = "white", lwd = 1.5)

<img src="README-quick_start-1.png" style="display: block; margin: auto;" />
