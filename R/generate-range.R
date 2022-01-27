#' Generate range polygons from an abundance raster
#'
#' Generate polygons encompassing the non-zero extent of a raster stack,
#' typically one giving the abundance of a species for different seasons or
#' weeks. In addition, polygons defining the extent of prediction area (i.e. the
#' region where the raster has some non-missing value, including zero) will be
#' generated.
#'
#' @param x `RasterStack`; the input raster to convert to polygons.
#' @param smooth logical; whether to smooth the polygons and drop small
#'   fragments consisting of single raster cells. Smoothing is done via the
#'   `smoothr` package. If `smooth = TRUE`, but the raw and smoothed polygons
#'   are returned as separate features.
#' @param drop_fill logical; whether to drop small polygon fragments and fill
#'   small holes in polygons.
#'
#' @details  The raw polygons will have sharp edges since they're being
#' converted from a square raster grid; however, by using `smooth = TRUE`, a set
#' of smooth polygons will be generated to more closely resemble the range
#' polygons seen in file guides. **Warning:** polygonizing the highest
#' resolution abundance data will take an extremely long time, it's best to use
#' the medium or low resolution data.
#'
#' @return An `sf` object containing the polygons defining the range. The range
#'   prediction area will each be encoded as a single feature distinguished by
#'   the `type` attribute. In addition, if `smooth = TRUE`, smoothed versions
#'   of these features will be generated and distinguished by the `smoothed`
#'   attribute. The `layer` attribute will contain the name of the raster
#'   layer that has been converted to polygons.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load seasonal abundance
#' abd <- load_raster(path, product = "abundance_seasonal", resolution = "lr")
#'
#' # convert to ranges
#' range_polygons <- generate_range(abd)
#' }
generate_range <- function(x, smooth = TRUE, drop_fill = TRUE) {
  if (is.character(x) || inherits(x, "SpatRaster")) {
    x <- raster::stack(x)
  }
  stopifnot(inherits(x, "Raster"))
  stopifnot(is.logical(smooth), length(smooth) == 1, !is.na(smooth))
  stopifnot(is.logical(drop_fill), length(drop_fill) == 1, !is.na(drop_fill))

  # suggests check
  if (!requireNamespace("smoothr", quietly = TRUE)) {
    stop("Package smoothr required for generate_range()")
  }
  if (!requireNamespace("units", quietly = TRUE)) {
    stop("Package units required for generate_range()")
  }

  # convert to polygons
  rng <- NULL
  pa <- NULL
  for (l in names(x)) {
    rng <- dplyr::bind_rows(rng, make_range_polygon(x[[l]],
                                                    type = "range"))
    pa <- dplyr::bind_rows(pa, make_range_polygon(x[[l]],
                                                  type = "prediction_area"))
  }
  rng$layer <- names(x)
  pa$layer <- names(x)
  rng$smoothed <- FALSE
  pa$smoothed <- FALSE

  # clean & smooth
  rng_pa <- rbind(rng, pa)
  if (isTRUE(smooth)) {
    # drop holes and polygons smaller than 1.5 times the cell size
    if (isTRUE(drop_fill)) {
      ca <- units::set_units(1.5 * prod(raster::res(x)), "m^2")
      rng_smooth <- smoothr::drop_crumbs(rng_pa, threshold = ca)
      rng_smooth <- smoothr::fill_holes(rng_smooth, threshold = ca)
    }

    # smooth
    rng_smooth <- smoothr::smooth(rng_smooth, method = "ksmooth",
                                  smoothness = 2, n = 5L)
    rng_smooth <- sf::st_make_valid(rng_smooth)

    # combine
    rng_smooth$smoothed <- TRUE
    rng_pa <- rbind(rng_pa, rng_smooth)
  }

  return(rng_pa)
}

make_range_polygon <- function(x, type = c("range", "prediction_area")) {
  stopifnot(inherits(x, "RasterLayer"), raster::nlayers(x) == 1)
  type <- match.arg(type)

  if (type == "range") {
    fun <- function(y) {
      y > 0
    }
  } else {
    fun <- function(y) {
      !is.na(y)
    }
  }

  p <- raster::rasterToPolygons(x, fun = fun, digits = 6)
  if (is.null(p)) {
    # empty polygon
    p <- sf::st_polygon()
    p <- sf::st_sfc(p, crs = sf::st_crs(x))
  } else {
    p <- p[!is.na(p@data[[1]]), ]
    p <- sf::st_as_sfc(p)
    p <- sf::st_set_precision(p, 1e4)
    p <- sf::st_union(p)
  }
  if (length(p) == 0) {
    # empty polygon
    p <- sf::st_polygon()
    p <- sf::st_sfc(p, crs = sf::st_crs(x))
  }
  p <- sf::st_sf(geometry = p)
  p$type <- type
  return(p)
}
