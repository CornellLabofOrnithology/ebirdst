#' Generate polygons for eBird Status and Trends stixels
#'
#' eBird Status and Trends divides space and time into variably sized "stixels"
#' within which individual base models are fit. The process of stixelization is
#' performed many times and the prediction at any given point is the median of
#' the predictions from all the stixels that that point falls in.
#' [load_stixels()] loads information on all the stixels that compromise a
#' species' Status and Trends model, with stixels identified by the location of
#' their centroid. This function uses this information to define polygons for
#' each stixel and attaches them to the original data in the form of an [sf]
#' object.
#'
#' @param x `data.frame` or [sf] object; stixel summary data loaded with
#'   [load_stixels()], or any other data frame with fields `lon`, `lat`,
#'   `stixel_width`, and `stixel_hight`.
#'
#' @return [sf] object with geometry column storing polygons representing the
#'   stixels boundaries.
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load stixel summary information
#' stixels <- load_stixels(path)
#'
#' # build stixel polygons
#' stixelize(stixels)
#' }
stixelize <- function(x) {
  UseMethod("stixelize")
}


#' @export
stixelize.data.frame <- function(x) {
  stopifnot(all(c("longitude", "latitude",
                  "stixel_width", "stixel_height") %in% names(x)))

  # drop empty stixels
  x <- x[x$stixel_width > 0 & x$stixel_height > 0, ]

  # function to make a single stixel
  f <- function(lon, lat, w, h) {
    bb <- sf::st_bbox(c(xmin = max(lon - w / 2, -180),
                        xmax = min(lon + w / 2, 180),
                        ymin = max(lat - h / 2, -90),
                        ymax = min(lat + h / 2, 90)))
    sf::st_as_sfc(bb)
  }
  stx <- sf::st_sfc(mapply(f,
                           x$longitude, x$latitude,
                           x$stixel_width, x$stixel_height),
                    crs = 4326)

  # combine with data
  sf::st_make_valid(sf::st_sf(x, geometry = stx))
}

#' @export
stixelize.sf <- function(x) {
  stopifnot(all(c("stixel_width", "stixel_height") %in% names(x)))
  stopifnot(all(sf::st_geometry_type(x) == "POINT"))

  x <- sf::st_transform(x, crs = 4326)

  # ensure coordinates are in df
  ll <- sf::st_coordinates(x)
  if (!all(c("longitude", "latitude") %in% names(x))) {
    x$longitude <- ll[, 1, drop = TRUE]
    x$latitude <- ll[, 2, drop = TRUE]
  }
  stixelize.data.frame(sf::st_set_geometry(x, NULL))
}
