#' Calculates spatial extent of non-zero data from Raster* object for plotting
#'
#' After loading a RasterStack of results, there are lots of NA values
#' and plots of individual raster layers will display at the full extent of the
#' study extent. To show an ideal extent, this function trims away 0 and NA
#' values and checks to make sure it returns a reasonable extent for plotting.
#' The returned Extent object can then be used for plotting. To access a
#' pre-calculated extent for the full annual cycle use
#' [load_fac_map_parameters()].
#'
#' @param x Raster* object; either full RasterStack or subset.
#'
#' @return raster Extent object
#'
#' @export
#'
#' @examples
#' \dontshow{
#' # simple toy example
#' r <- raster::raster(nrow = 100, ncol = 100)
#' r[5025:5075] <- 1
#' raster::extent(r)
#' calc_full_extent(r)
#' }
#' \donttest{
#' # download and load example abundance data
#' sp_path <- ebirdst_download("example_data")
#' abd <- load_raster("abundance", sp_path)
#'
#' # calculate full extent
#' plot_extent <- calc_full_extent(abd)
#'
#' # plot
#' raster::plot(abd[[1]], axes = FALSE, ext = plot_extent)
#' }
calc_full_extent <- function(x) {
  stopifnot(inherits(x, "Raster"))

  # aggregate stack for speed, otherwise everything else takes too long
  stack <- raster::aggregate(x, fact = 3)

  # convert 0s to NAs, otherwise trimming is slow and the extent is too broad
  stack[stack == 0] <- NA

  # work on individual layers because trim has a bug when applied to stacks
  # also this approach is faster
  e <- c(NA_real_, NA_real_, NA_real_, NA_real_)
  for (i in seq_len(raster::nlayers(stack))) {
    e_trim <- raster::extent(raster::trim(stack[[i]]))
    e[1] <- min(e[1], e_trim[1], na.rm = TRUE)
    e[2] <- max(e[2], e_trim[2], na.rm = TRUE)
    e[3] <- min(e[3], e_trim[3], na.rm = TRUE)
    e[4] <- max(e[4], e_trim[4], na.rm = TRUE)
  }
  e <- raster::extent(e)

  # sometimes extent calculations get weird and you'll get a very broad
  # extent that goes further than you want, so check against the input
  e_input <- raster::extent(x)
  e[1] <- max(e[1], e_input[1], na.rm = TRUE)
  e[2] <- min(e[2], e_input[2], na.rm = TRUE)
  e[3] <- max(e[3], e_input[3], na.rm = TRUE)
  e[4] <- min(e[4], e_input[4], na.rm = TRUE)

  return(e)
}


#' Calculates bins (breaks) based on standard deviations of Box-Cox
#' power-transformed data for mapping
#'
#' Mapping species abundance across the full-annual cycle presents a challenge,
#' in that patterns of concentration and dispersion in abundance change
#' throughout the year, making it difficult to define color bins that suit all
#' seasons and accurately reflect the detail of abundance predictions. To
#' address this, we selected a method (described by Maciejewski et al. 2013)
#' that first selects an optimal power (the Box-Cox method) for normalizing
#' the data, then power transforms the entire year of non-zero data, constructs
#' bins with the power-transformed data using standard-deviations, and then
#' un-transforms the bins. To access a pre-calculated bins for the full annual
#' cycle use [load_fac_map_parameters()].
#'
#' @param x RasterStack or RasterBrick; original eBird Status and Trends product
#'   raster GeoTIFF with 52 bands, one for each week.
#'
#' @return A list with two elements: `bins` is a vector containing the break
#'   points of the bins and `power` is the optimal power used to transform data
#'   when calculating bins.
#'
#' @export
#'
#' @references Ross Maciejewski, Avin Pattah, Sungahn Ko, Ryan Hafen, William S.
#' Cleveland, David S. Ebert.  Automated Box-Cox Transformations for Improved
#' Visual Encoding. IEEE Transactions on Visualization and Computer Graphics,
#' 19(1): 130-140, 2013.
#'
#' @examples
#' # download and load example abundance data
#' sp_path <- ebirdst_download("example_data")
#' abd <- load_raster("abundance", sp_path)
#' \dontshow{
#' # crop to speed up cran tests
#' e <-  raster::extent(abd)
#' e[2] <- e[1] + (e[2] - e[1]) / 4
#' e[4] <- e[3] + (e[4] - e[3]) / 4
#' abd <- raster::crop(abd[[30]], e)
#' }
#'
#' # calculate bins for a single week for this example
#' year_bins <- calc_bins(abd)
calc_bins <- function(x) {
  stopifnot(inherits(x, "Raster"))
  if (all(is.na(suppressWarnings(raster::maxValue(x)))) &
      all(is.na(suppressWarnings(raster::minValue(x))))) {
    stop("Input Raster* object must have non-NA values.")
  }
  if (all(raster::maxValue(x) == 0)) {
    stop("Raster must have at least 2 non-zero values to calculate bins")
  }

  # get a vector of all the values in the stack
  zrv <- raster::getValues(x)
  vals_for_pt <- zrv[!is.na(zrv) & zrv > 0]
  if (length(vals_for_pt) <= 1) {
    stop("Raster must have at least 2 non-zero values to calculate bins")
  }

  # box-cox transform
  pt <- car::powerTransform(vals_for_pt)
  this_power <- pt$lambda

  lzwk <- vals_for_pt ^ this_power
  rm(zrv)

  # setup the binning structure
  # calculate metrics
  maxl <- max(lzwk, na.rm = TRUE)
  minl <- min(lzwk, na.rm = TRUE)
  mdl <- mean(lzwk, na.rm = TRUE)
  sdl <- stats::sd(lzwk, na.rm = TRUE)
  rm(lzwk)

  # build a vector of bins
  log_sd <- c(mdl - (3.00 * sdl), mdl - (2.50 * sdl), mdl - (2.00 * sdl),
              mdl - (1.75 * sdl), mdl - (1.50 * sdl), mdl - (1.25 * sdl),
              mdl - (1.00 * sdl), mdl - (0.75 * sdl), mdl - (0.50 * sdl),
              mdl - (0.25 * sdl), mdl - (0.125 * sdl),
              mdl,
              mdl + (0.125 * sdl), mdl + (0.25 * sdl),
              mdl + (0.50 * sdl), mdl + (0.75 * sdl), mdl + (1.00 * sdl),
              mdl + (1.25 * sdl), mdl + (1.50 * sdl), mdl + (1.75 * sdl),
              mdl + (2.00 * sdl), mdl + (2.50 * sdl), mdl + (3.00 * sdl))

  # lots of checks for values outside of the upper and lower bounds

  # remove +3 sd break if it is greater than max
  if (maxl < mdl + (3.00 * sdl)) {
    log_sd <- log_sd[1:length(log_sd) - 1]
  }

  # add max if the max is greater than +3 sd break
  if (maxl > mdl + (3.00 * sdl) | maxl > log_sd[length(log_sd)]) {
    log_sd <- append(log_sd, maxl)
  }

  # remove the -3 sd break if it is less than the min
  if (minl > mdl - (3.00 * sdl)) {
    log_sd <- log_sd[2:length(log_sd)]
  }

  # add min if the min is less than -3 sd break
  if (minl < mdl - (3.00 * sdl) | minl < log_sd[1]) {
    log_sd <- append(log_sd, minl, after = 0)
  }

  if (log_sd[1] < 0) {
    log_sd[1] <- 0.01 ^ this_power
  }

  if (log_sd[1] ^ (1 / this_power) < 0.01) {
    log_sd[1] <- 0.01 ^ this_power
  }

  # untransform
  bins <- log_sd ^ (1 / this_power)
  rm(log_sd)

  # if transform power was negative, flip bins
  if (this_power < 0) {
    bins <- rev(bins)
  }

  return(list(bins = sort(unname(bins)), power = unname(this_power)))
}


#' Map PI centroid locations
#'
#' Creates a map showing the stixel centroid locations for predictor importance
#' values, with an optional spatiotemporal subset using an [ebirdst_extent]
#' object
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal extent to
#'   filter the data to.
#'
#' @return Plot showing locations of PI centroids.
#'
#' @export
#'
#' @examples
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # define a spatiotemporal extent to plot
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' map_centroids(path = sp_path, ext = e)
map_centroids <- function(path, ext) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  if (missing(ext)) {
    stop("A spatiotemporal extent must be provided.")
  } else {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # colors
  pal <- list(pd_a = "#1b9377", pd_s = "#b3e2cd",
              pi_a = "#d95f02", pi_s = "#fdcdac")
  pal <- list(pd_a = "#31a354", pd_s = "#a1d99b",
              pi_a = "#d95f0e", pi_s = "#feb24c")

  # convert data to spatial, find bounding box
  # pis
  pis <- load_pis(path = path)
  pis <- dplyr::distinct(pis[, c("lon", "lat", "date")])
  pis <- sf::st_as_sf(pis, coords = c("lon", "lat"), crs = 4326)
  pis <- sf::st_transform(pis, crs = mollweide)
  # bbox
  bb <- sf::st_as_sfc(sf::st_bbox(pis))

  # initialize graphical parameters
  p <- graphics::par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), bg = "white")

  # plot base map
  graphics::plot(bb, col = NA, border = NA)
  graphics::plot(sf::st_geometry(ned_wh_co_moll),
                 col = "#5a5a5a", border = "#222222", lwd = 1, add = TRUE)
  # label setup
  usr <- graphics::par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]

  # plotting pis
  # first plot all possible pis
  graphics::plot(sf::st_geometry(pis), col = pal$pi_a, cex = 0.4, pch = 16,
                 add = TRUE)
  graphics::text(x = usr[1] + 0.05 * xwidth,
                 y = usr[3] + 0.12 * yheight,
                 adj = 0,
                 paste("Available PIs: ", nrow(pis), sep = ""),
                 cex = 1,
                 col = pal$pi_a)

  # plot pis within extent
  if (!missing(ext)) {
    pis_sub <- ebirdst_subset(pis, ext)
    graphics::plot(sf::st_geometry(pis_sub),
                   col = pal$pi_s, cex = 0.4, pch = 16,
                   add = TRUE)
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.08 * yheight,
                   adj = 0,
                   paste("Selected PIs: ", nrow(pis_sub), sep = ""),
                   cex = 1,
                   col = pal$pi_s)
  }

  # plot reference data
  graphics::plot(sf::st_geometry(ned_wh_co_moll),
                 col = NA, border = "#222222", lwd = 1, add = TRUE)
  graphics::plot(sf::st_geometry(ned_wh_st_moll),
                 col = NA, border = "#222222", lwd = 0.75, add = TRUE)

  # label
  usr <- graphics::par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]

  # pis
  graphics::text(x = usr[1] + 0.05 * xwidth,
                 y = usr[3] + 0.12 * yheight,
                 adj = 0,
                 paste("Available PIs: ", nrow(pis), sep = ""),
                 cex = 1,
                 col = pal$pi_a)
  if (!missing(ext)) {
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.08 * yheight,
                   adj = 0,
                   paste("Selected PIs: ", nrow(pis_sub), sep = ""),
                   cex = 1,
                   col = pal$pi_s)
  }

  graphics::par(p)
  invisible()
}


#' Calculate and map effective extent of selected centroids
#'
#' The selection of stixel centroids for analysis of predictor importances (PIs) yields an
#' effective footprint, or extent, showing the effective location of where the
#' information going into the analysis with PIs is based. While a
#' bounding box or polygon may be used to select a set of centroids, due to the
#' models being fit within large rectangular areas, the information from a set
#' of centroids often comes from the core of the selected area. This function
#' calculates where the highest proportion of information is coming from,
#' returns a raster and plots that raster, with the selected area and centroids
#' for reference. The legend shows, for each pixel, what percentage of the
#' selected stixels are contributing information, ranging from 0 to 1.
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal
#'   extent to filter the data to.
#' @param plot logical; whether to plot the results or just return a raster
#'   without plotting.
#'
#' @return A raster showing the percentage of the selected stixels that are
#'   contributing to each grid cell. In addition, if `plot = TRUE` this raster
#'   will be plotted along with centroid locations and [ebirdst_extent]
#'   boundaries.
#'
#' @export
#'
#' @examples
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # define a spatioremporal extent
#' bb_vec <- c(xmin = -86, xmax = -84, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' \donttest{
#' # calculate effective extent map
#' eff <- calc_effective_extent(path = sp_path, ext = e)
#' }
calc_effective_extent <- function(path, ext, plot = TRUE) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(inherits(ext, "ebirdst_extent"))
  stopifnot(is.logical(plot), length(plot) == 1)
  if (all(c(0, 1) == round(ext$t, 2))) {
    warning(paste("Without temporal limits in ext, this function will take",
                  "considerably longer to run and is less informative."))
  }

  # load data
  r_tmplt <- load_raster(product = "template", path = path)
  pipd <- load_pis(path = path)

  # subset
  pipd <- dplyr::distinct(pipd[, c("lon", "lat", "date", "stixel_width",
                                   "stixel_height")])
  pipd <- ebirdst_subset(pipd, ext = ext)
  # stixelize
  stixels <- stixelize(pipd)
  # project to template raster projection
  stixels <- sf::st_transform(stixels, crs = sf::st_crs(r_tmplt))

  # summarize: % of stixels overlapping each cell
  r_stix <- raster::crop(
    fasterize::fasterize(stixels, r_tmplt, fun = "count"),
    raster::extent(stixels))
  r_stix <- r_stix / nrow(stixels)

  # plot
  if (isTRUE(plot)) {
    r_stix_moll <- suppressWarnings(raster::projectRaster(r_stix,
                                                          crs = mollweide,
                                                          method = "ngb"))
    r_stix_moll[r_stix_moll >= 1] <- 1
    stixels_moll <- sf::st_transform(stixels, crs = mollweide)
    pipd_sf <- sf::st_as_sf(pipd, coords = c("lon", "lat"), crs = 4326)
    pipd_sf <- sf::st_transform(pipd_sf, mollweide)

    # convert extent to polygon and mollweide for plotting
    if (ext$type == "bbox") {
      ext_poly <- sf::st_as_sfc(ext$extent)
    } else if (ext$type == "polygon") {
      ext_poly <- ext$extent
    } else {
      stop("Spatiotemporal extent type not accepted.")
    }
    ext_poly_moll <- sf::st_transform(ext_poly, crs = mollweide)

    # plot
    p <- graphics::par(mar = c(0.25, 0.25, 0.25, 0.25), bg = "#ffffff")

    raster::plot(r_stix_moll, ext = raster::extent(stixels_moll),
                 breaks = c(0, seq(0.5, 1, by = 0.05)),
                 col = viridisLite::viridis(11), colNA = "#000000",
                 maxpixels = raster::ncell(r_stix_moll),
                 axis.args = list(at = c(0, seq(0.5, 1, by = 0.1)),
                                  labels = c(0, seq(0.5, 1, by = 0.1))),
                 axes = FALSE, box = FALSE, legend = TRUE)

    graphics::plot(sf::st_geometry(ned_wh_co_moll),
                   border = "#ffffff", col = NA, lwd = 1.5, add = TRUE)
    graphics::plot(sf::st_geometry(ned_wh_st_moll),
                   border = "#ffffff", col = NA, lwd = 1, add = TRUE)
    graphics::plot(sf::st_geometry(ext_poly_moll),
                   border = "red", col = NA, lwd = 1.5, add = TRUE)
    graphics::plot(sf::st_geometry(pipd_sf),
                   col = "#000000", pch = 16, cex = 1 * graphics::par()$cex,
                   add = TRUE)
    graphics::par(p)
  }
  invisible(r_stix)
}


#' eBird Status and Trends color palettes for abundance data
#'
#' Generate the color palettes used in the eBird Status and Trends relative
#' abundance maps.
#'
#' @param n integer; the number of colors to be in the palette.
#' @param season character; the season to generate colors for or "weekly" to
#'   get the color palette used in the weekly abundance animations.
#'
#' @return A character vector of color hex codes.
#' @export
#'
#' @examples
#' # breeding season color palette
#' abundance_palette(10, season = "breeding")
abundance_palette <- function(n, season = c("weekly",
                                            "breeding", "nonbreeding",
                                            "migration",
                                            "prebreeding_migration",
                                            "postbreeding_migration",
                                            "year_round")) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  season <- match.arg(season)

  # set base color by season
  col_zero <- "#dddddd"
  if (season == "weekly") {
    plsm <- rev(viridisLite::plasma(n - 1, end = 0.9))
    gry <- grDevices::colorRampPalette(c(col_zero, plsm[1]))
    return(c(gry(4)[2], plsm))
  } else if (season == "breeding") {
    base_col <- "#cc503e"
  } else if (season == "nonbreeding") {
    base_col <- "#1d6996"
  } else if (season %in% c("migration", "postbreeding_migration")) {
    base_col <- "#edad08"
  } else if (season == "prebreeding_migration") {
    base_col <- "#73af48"
  } else if (season == "year_round") {
    base_col <- "#6f4070"
  } else {
    stop("Invalid season.")
  }

  # seasonal palettes
  gry <- grDevices::colorRampPalette(c(col_zero, base_col))
  mid <- grDevices::colorRampPalette(c(gry(5)[2], base_col))
  black <- grDevices::colorRampPalette(c(base_col, "#000000"))
  pal <- grDevices::colorRampPalette(c(gry(5)[2], mid(9)[5], base_col,
                                       black(5)[2]))
  return(pal(n))
}
