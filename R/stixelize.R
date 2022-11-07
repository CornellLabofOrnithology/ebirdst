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
#' @param x data frame; stixel summary data loaded with [load_stixels()], or any
#'   other data frame with fields `lonitude_min`, `lontidue_max`,
#'   `latitude_min`, and `latitude_max`.
#'
#' @return [sf] object with geometry column storing polygons representing the
#'   stixels boundaries.
#' @export
#'
#' @examples
#' \dontrun{
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
  stopifnot(all(c("longitude_min", "longitude_max",
                  "latitude_min", "latitude_max") %in% names(x)))

  # function to make a single stixel
  f <- function(lon_min, lon_max, lat_min, lat_max) {
    # if stixel crosses dateline, split it into two
    if (lon_min < 0 && lon_max > 0) {
      bb_west <- sf::st_bbox(c(xmin = -179.9999, xmax = lon_min,
                               ymin = lat_min, ymax = lat_max))
      bb_east <- sf::st_bbox(c(xmin = lon_max, xmax = 179.9999,
                               ymin = lat_min, ymax = lat_max))
      bb <- c(sf::st_as_sfc(bb_west), sf::st_as_sfc(bb_east))
      bb <- sf::st_combine(bb)
    } else {
      bb <- sf::st_bbox(c(xmin = lon_min, xmax = lon_max,
                          ymin = lat_min, ymax = lat_max))
      bb <- sf::st_as_sfc(bb)
    }
    return(sf::st_make_valid(bb))
  }
  stx <- sf::st_sfc(mapply(f,
                           x$longitude_min, x$longitude_max,
                           x$latitude_min, x$latitude_max),
                    crs = 4326)

  # combine with data
  sf::st_sf(x, geometry = stx)
}
