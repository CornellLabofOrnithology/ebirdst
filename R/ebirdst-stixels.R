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
  stopifnot(all(c("lon", "lat", "stixel_width", "stixel_height") %in% names(x)))

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
  stx <- sf::st_sfc(mapply(f, x$lon, x$lat, x$stixel_width, x$stixel_height),
                    crs = 4326)

  # combine with data
  sf::st_sf(x, geometry = stx)
}

#' @export
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


#' Calculate the spatial footprint of a set of stixels
#'
#' eBird Status and Trends divides space and time into variably sized "stixels"
#' within which individual base models are fit. The process of stixelization is
#' performed many times and the prediction at any given point is the median of
#' the predictions from all the stixels that that point falls in. For a given
#' spatiotemporal extent, this function identifies the set of stixels whose
#' centroids fall within that extent and calculates the spatial footprint of
#' these stixels, i.e. a surface indicating the proportion of the selected
#' stixels that contribute information to model estimates at each location. This
#' footprint gives an estimate of where the information for the model
#' predictions, predictor importances (PIs), and partial dependencies (PDs) come
#' from.
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to filter the
#'   data to.
#'
#' @return A [stixel_footprint] object consisting of a list with three elements:
#'   - `footprint`: a `RasterStack` giving the percentage of the selected
#'   stixels that are contributing to each grid cell.
#'   - `centroids`: an [sf] object containing the stixel centroids points.
#'   - `extent`: an [ebirdst_extent] object specifying the chosen spatiotemporal
#'   extent.
#'
#'   The stixel footprint can be mapped by calling [plot()] on the returned
#'   [stixel_footprint] object.
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
#' # define a spatiotemporal extent
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 42, ymax = 45)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' # calculate effective extent map
#' footprint <- stixel_footprint(path, ext = e)
#' plot(footprint)
#' }
stixel_footprint <- function(path, ext) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(inherits(ext, "ebirdst_extent"))
  if (all(c(0, 1) == round(ext$t, 2))) {
    warning(paste("Without temporal limits in ext, this function will take",
                  "considerably longer to run and is less informative."))
  }

  # load data
  r_tmplt <- load_raster(product = "template", path = path)
  stx <- load_stixels(path = path)

  # subset
  stx <- ebirdst_subset(stx, ext = ext)
  centroids <- sf::st_as_sf(stx, coords = c("lon", "lat"), crs = 4326)
  # stixelize
  stx_sf <- stixelize(stx)
  # project to template raster projection
  stx_sf <- sf::st_transform(stx_sf, crs = sf::st_crs(r_tmplt))

  # summarize: % of stixels overlapping each cell
  stx_r <- fasterize::fasterize(stx_sf, r_tmplt, fun = "sum")
  stx_r <- raster::crop(stx_r, raster::extent(stx_sf))
  stx_r <- stx_r / nrow(stx_sf)
  stx_r[stx_r >= 1] <- 1

  structure(list(footprint = stx_r, centroids = centroids, extent = ext),
            class = "stixel_footprint")
}

#' @param x [stixel_footprint] object to map.
#' @param ... ignored.
#' @export
#' @rdname stixel_footprint
plot.stixel_footprint <- function(x, ...) {
  # convert extent to polygon
  if (x$extent$type == "bbox") {
    ext_poly <- sf::st_as_sfc(x$extent$extent)
  } else if (x$extent$type == "polygon") {
    ext_poly <- x$extent$extent
  } else {
    stop("Spatiotemporal extent type not accepted.")
  }

  # project to eck4
  stx_r_eck <- suppressWarnings(raster::projectRaster(x$footprint,
                                                      crs = prj_eck4,
                                                      method = "ngb"))
  centroids <- sf::st_transform(x$centroids, prj_eck4)
  ext_poly_eck <- sf::st_transform(ext_poly, crs = prj_eck4)
  # plot
  p <- graphics::par(mar = c(0.25, 0.25, 0.25, 0.25), bg = "#ffffff")

  raster::plot(stx_r_eck, ext = calc_full_extent(stx_r_eck),
               breaks = seq(0, 1, by = 0.05),
               col = viridisLite::viridis(20), colNA = "#000000",
               maxpixels = raster::ncell(stx_r_eck),
               axis.args = list(at = seq(0, 1, by = 0.25),
                                labels = seq(0, 1, by = 0.25)),
               axes = FALSE, box = FALSE, legend = TRUE)

  graphics::plot(sf::st_geometry(ne_adm0_eck),
                 border = "#ffffff", col = NA, lwd = 1.5, add = TRUE)
  graphics::plot(sf::st_geometry(ne_adm1_eck),
                 border = "#ffffff", col = NA, lwd = 1, add = TRUE)
  graphics::plot(sf::st_geometry(ext_poly_eck),
                 border = "red", col = NA, lwd = 1.5, add = TRUE)
  graphics::plot(sf::st_geometry(centroids),
                 col = "#000000", pch = 16, cex = 1 * graphics::par()$cex,
                 add = TRUE)
  graphics::par(p)
  invisible(NULL)
}
