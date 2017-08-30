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
#' str(template_raster)
#'
#' r <- raster()
#' extended_raster <- extend(r, extent(template_raster))
"template_raster"
