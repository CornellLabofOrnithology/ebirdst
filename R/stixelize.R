#' Generate stixel polygons from PI data
#'
#' All predictor importance data are provided at the stixel level. In these
#' files, the stixel is defined based on a centroid, width, and height. This
#' function uses this information to define polygons for each stixel and
#' attaches them to the original data in the form of an [sf] object
#'
#' @param x `data.frame` or [sf] object; PI data loaded with [load_pis()], or
#'   any other data frame with fields `lon`, `lat`, `stixel_width`, and
#'   `stixel_hight`.
#'
#' @return [sf] object with geometry column storing polygons representing the s
#'   tixels boundaries.
#' @export
#'
#' @examples
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # load predictor importance
#' pis <- load_pis(sp_path)
#'
#' stixelize(pis)
#'
#' # also works on sf objects
#' pis_sf <- sf::st_as_sf(pis, coords = c("lon", "lat"), crs = 4326)
#' stixelize(pis_sf)
stixelize <- function(x) {
  UseMethod("stixelize")
}


#' @export
#' @describeIn stixelize PI or PD data
stixelize.data.frame <- function(x) {
  stopifnot(all(c("lon", "lat", "stixel_width", "stixel_height") %in% names(x)))

  # function to make a single stixel
  f <- function(lon, lat, w, h) {
    bb <- sf::st_bbox(c(xmin = lon - w / 2, xmax = lon + w / 2,
                        ymin = lat - h / 2, ymax = lat + h / 2))
    sf::st_as_sfc(bb)
  }
  stx <- sf::st_sfc(mapply(f, x$lon, x$lat, x$stixel_width, x$stixel_height),
                    crs = 4326)

  # combine with data
  sf::st_sf(x, geometry = stx)
}


#' @export
#' @describeIn stixelize PI or PD data as `sf` object
stixelize.sf <- function(x) {
  stopifnot(all(c("stixel_width", "stixel_height") %in% names(x)))
  stopifnot(all(sf::st_geometry_type(x) == "POINT"))

  x <- sf::st_transform(x, crs = 4326)

  # ensure coordinates are in df
  ll <- sf::st_coordinates(x)
  if (!all(c("lon", "lat") %in% names(x))) {
    x$lon <- ll[, 1, drop = TRUE]
    x$lat <- ll[, 2, drop = TRUE]
  }
  stixelize.data.frame(sf::st_set_geometry(x, NULL))
}
