#' Template raster for extending the extents of STEM raster layers
#'
#' A raster layer that is used to unify the extents of STEM results across weeks
#' of the year. Currently, each individual STEM result layer (e.g.,
#' "abundance_umean") is stored by having all NA values trimmed, which results
#' in raster layers that have different extents for each species for each week
#' of the year. This template raster provides the full spatial extent of the
#' study area.
#'
#' @format A RasterLayer with pixels having value of 0 that represents the full
#' spatial extent, the resolution, and the coordinate system of the STEM study
#' extent.
#'
#' @examples
#' str(stemhelper::template_raster)
#'
#' r <- raster()
#' extended_raster <- extend(r, extent(stemhelper::template_raster))
"template_raster"

#' RasterStack of zero values based on checklist-level ensemble support
#'
#' A RasterStack (one layer for each week of the year) indicating the ensemble support based on sufficient checklists to fit models within an individiual stixel. Values range from 0 to 100. For presentation, we use a value of 99 to indicate sufficient information to assume that, if there are no positive occurrence or abundance predictions within a given raster cell, that there is a value of 0 occurrence and 0 abundance.
#'
#' @format A RasterStack, where each layer is a week of the year and where the values range of 0 to 100 representing the number of contributing independent models (stixels) with sufficient checklists to fit and predict within that indepdendent model (stixel). Each RasterLayer in the RasterStack has the full spatial extent and same resolution and coordinate system of the STEM study extent.
#'
#' @examples
#' stemhelper::zero_es_stack
#'
#' # create binary layers for each week where to show assumed zeroes
#' stemhelper::zero_es_stack[stemhelper::zero_es_stack < 99] <- NA
"zero_es_stack"
