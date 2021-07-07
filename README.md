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

Install `ebirdst` from CRAN with:

    install.packages("ebirdst")

Alternatively, you can install the development version from GitHub with:

    # install.packages("remotes")
    remotes::install_github("CornellLabofOrnithology/ebirdst")

The current version of `ebirdst` is designed for working with the eBird
Status and Trends data products released in December 2020, with
abundance estimates made for 2019. **Users are strongly discouraged from
comparing Status and Trends results between years.**

## Citation

If you use the the eBird Status & Trends data please cite it with:

<blockquote>
Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S.
Ligocki, W. Hochachka, C. Wood, I. Davies, M. Iliff, L. Seitz. 2020.
eBird Status and Trends, Data Version: 2019; Released: 2020, Cornell Lab
of Ornithology, Ithaca, New York.
<a href="https://doi.org/10.2173/ebirdst.2019" class="uri">https://doi.org/10.2173/ebirdst.2019</a>
</blockquote>

## Vignettes

For full package documentation, including a series of vignettes covering
the full spectrum from introductory to advanced usage, please see the
package [website](https://cornelllabofornithology.github.io/ebirdst).
There are four available vignettes are:

-   [Introduction](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-introduction.html):
    details data access and the structure of the results.
-   [Introductary
    mapping](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-intro-mapping.html):
    covers the basic of making abundance maps using the eBird Status and
    Trends data.
-   [Advanced
    mapping](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-advanced-mapping.html):
    shows how to reproduce the seasonal abundance maps, range polygons,
    and statistics on the [eBird Status and Trends
    website](https://ebird.org/science/status-and-trends).
-   [Non-raster
    data](https://cornelllabofornithology.github.io/ebirdst/articles/ebirdst-pipd.html):
    covers working with the non-raster data in the `ebirst` data
    packages, including exploring **predictor importance**, plotting
    **partial dependence** curves for model predictors, and assessing
    model quality with **predictive performance metrics**.

To access old versions of the vignettes suitable for working with
previous eBird Status and Trends releases go to the [versions page of
the GitHub
repository](https://github.com/CornellLabofOrnithology/ebirdst/releases).
Download and unzip the source code for the release associated with the
version of the data you want to work with. Opening the file
`docs/index.html` will open the package website for this version, which
gives access to all the documentation and vignettes.

## Quick Start

This quick start guide shows how to download data and plot abundance
values similar to how they are plotted for the [eBird Status and Trends
Abundance animated
maps](https://ebird.org/science/status-and-trends/yebsap/abundance-map-weekly).
In this guide, and throughout all package documentation, a simplified
example dataset is used consisting of Yellow-bellied Sapsucker in
Michigan. For a full list of the species available for download, look at
the data frame `ebirst_runs`, which is included in this package.

**Important note: after downloading the results, do not change the file
structure.** All functionality in this package relies on the structure
inherent in the delivered results. Changing the folder and file
structure will cause errors with this package.

    library(ebirdst)
    library(raster)
    library(sf)
    library(fields)
    library(rnaturalearth)

    # download example data, yellow-bellied sapsucker in michigan
    dl_path <- ebirdst_download(species = "example_data")

    # load relative abundance raster stack with 52 layers, one for each week
    abd <- load_raster("abundance", path = dl_path)

    # load species specific mapping parameters
    pars <- load_fac_map_parameters(dl_path)
    # custom coordinate reference system
    crs <- pars$custom_projection
    # legend breaks
    breaks <- pars$abundance_bins
    # legend labels for top, middle, and bottom
    labels <- attr(breaks, "labels")

    # get a date vector specifying which week each raster layer corresponds to
    weeks <- parse_raster_dates(abd)
    print(weeks)
    #>  [1] "2019-01-04" "2019-01-11" "2019-01-18" "2019-01-25" "2019-02-01" "2019-02-08" "2019-02-15" "2019-02-22" "2019-03-01"
    #> [10] "2019-03-08" "2019-03-15" "2019-03-22" "2019-03-29" "2019-04-05" "2019-04-12" "2019-04-19" "2019-04-26" "2019-05-03"
    #> [19] "2019-05-10" "2019-05-17" "2019-05-24" "2019-05-31" "2019-06-07" "2019-06-14" "2019-06-21" "2019-06-28" "2019-07-06"
    #> [28] "2019-07-13" "2019-07-20" "2019-07-27" "2019-08-03" "2019-08-10" "2019-08-17" "2019-08-24" "2019-08-31" "2019-09-07"
    #> [37] "2019-09-14" "2019-09-21" "2019-09-28" "2019-10-05" "2019-10-12" "2019-10-19" "2019-10-26" "2019-11-02" "2019-11-09"
    #> [46] "2019-11-16" "2019-11-23" "2019-11-30" "2019-12-07" "2019-12-14" "2019-12-21" "2019-12-28"

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

<img src="README-quick_start-1.png" style="display: block; margin: auto;" />
