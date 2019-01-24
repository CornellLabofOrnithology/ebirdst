#' Calculates spatial extent of non-zero data from Raster* object for plotting
#'
#' After loading a RasterStack of results, there are lots of NA values
#' and plots of individual raster layers will display at the full extent of the
#' study extent. To show an ideal extent, this function trims away 0 and
#' NA values and checks to make sure it returns a reasonable extent for
#' plotting. The returned Extent object can then be used for plotting.
#'
#' @param x Raster* object; either full RasterStack or subset.
#'
#' @return raster Extent object
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example abundance data
#' sp_path <- download_data("example_data")
#' abd <- load_raster("abundance_umean", sp_path)
#'
#' # calculate full extent
#' plot_extent <- calc_full_extent(abd)
#'
#' # plot
#' raster::plot(abd[[1]], axes = FALSE, ext = plot_extent)
#'
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
#' un-transforms the bins.
#'
#' @param x RasterStack or RasterBrick; original eBird Status and Trends product
#'   raster GeoTiff with 52 bands, one for each week.
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
#' \dontrun{
#'
#' # download and load example abundance data
#' sp_path <- download_data("example_data")
#' abd <- load_raster("abundance_umean", sp_path)
#'
#' # calculate bins
#' year_bins <- calc_bins(abd)
#'
#' # plot
#' raster::plot(abd[[30]], axes = FALSE, breaks = year_bins$bins,
#'              col = viridisLite::viridis(length(year_bins$bins) - 1))
#'
#' }
calc_bins <- function(x) {
  stopifnot(inherits(x, "Raster"))

  if (all(is.na(suppressWarnings(raster::maxValue(x)))) &
      all(is.na(suppressWarnings(raster::minValue(x))))) {
    stop("Input Raster* object must have non-NA values.")
  }

  # get a vector of all the values in the stack
  zrv <- raster::getValues(x)
  vals_for_pt <- zrv[!is.na(zrv) & zrv > 0]

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

  return(list(bins = bins, power = this_power))
}


#' Map PI and PD centroid locations
#'
#' Creates a map showing the stixel centroid locations for PIs and/or PDs, with
#' an optional spatiotemporal subset using an [ebirdst_extent] object
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal extent to
#'   filter the data to.
#' @param plot_pis logical; whether to show PI stixel centroid locations.
#' @param plot_pds logical; whether to show PD stixel centroid locations.
#'
#' @return Plot showing locations of PIs and/or PDs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example data
#' sp_path <- download_data("example_data")
#'
#' # define a spatiotemporal extent to plot
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' map_centroids(path = sp_path, ext = e)
#'
#' }
map_centroids <- function(path, ext, plot_pis = TRUE, plot_pds = TRUE) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(is.logical(plot_pis), length(plot_pis) == 1)
  stopifnot(is.logical(plot_pds), length(plot_pds) == 1)
  if (!plot_pds & !plot_pis) {
    stop("Plotting of both PIs and PDs set to FALSE. Nothing to plot!")
  }
  if (missing(ext)) {
    stop("A spatiotemporal extent must be provided.")
  } else {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # convert data to spatial, find bounding box
  # pds
  if (isTRUE(plot_pds)) {
    pds <- load_pds(path = path)
    pds <- dplyr::distinct(pds[, c("lon", "lat", "date")])
    pds <- sf::st_as_sf(pds, coords = c("lon", "lat"), crs = 4326)
    pds <- sf::st_transform(pds, crs = mollweide)
  } else {
    pds <- NULL
  }
  # pis
  if (isTRUE(plot_pis)) {
    pis <- load_pis(path = path)
    pis <- dplyr::distinct(pis[, c("lon", "lat", "date")])
    pis <- sf::st_as_sf(pis, coords = c("lon", "lat"), crs = 4326)
    pis <- sf::st_transform(pis, crs = mollweide)
  } else {
    pis <- NULL
  }
  # bbox
  bb <- sf::st_as_sfc(sf::st_bbox(rbind(pis, pds)))

  # initialize graphical parameters
  p <- graphics::par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), bg = "white")

  # plot base map
  graphics::plot(bb, col = NA, border = NA)
  graphics::plot(sf::st_geometry(ned_wh_co_moll),
                 col = "#5a5a5a", border = "#000000", lwd = 1, add = TRUE)
  # label setup
  usr <- graphics::par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]

  # plotting pds
  if (isTRUE(plot_pds)) {
    # first plot all possible pds
    graphics::plot(sf::st_geometry(pds), col = "#1b9377", cex = 0.4, pch = 16,
                   add = TRUE)
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.2 * yheight,
                   adj = 0,
                   paste("Available PDs: ", nrow(pds), sep = ""),
                   cex = 1,
                   col = "#1b9377")

    # plot pds within extent
    if (!missing(ext)) {
      pds_sub <- ebirdst_subset(pds, ext)
      graphics::plot(sf::st_geometry(pds_sub),
                     col = "#b3e2cd", cex = 0.4, pch = 16,
                     add = TRUE)
      graphics::text(x = usr[1] + 0.05 * xwidth,
                     y = usr[3] + 0.16 * yheight,
                     adj = 0,
                     paste("Selected PDs: ", nrow(pds_sub), sep = ""),
                     cex = 1,
                     col = "#b3e2cd")
    }
  }

  # plotting pis
  if (isTRUE(plot_pis)) {
    # first plot all possible pis
    graphics::plot(sf::st_geometry(pis), col = "#d95f02", cex = 0.4, pch = 16,
                   add = TRUE)
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.12 * yheight,
                   adj = 0,
                   paste("Available PIs: ", nrow(pis), sep = ""),
                   cex = 1,
                   col = "#d95f02")

    # plot pds within extent
    if (!missing(ext)) {
      pis_sub <- ebirdst_subset(pis, ext)
      graphics::plot(sf::st_geometry(pis_sub),
                     col = "#fdcdac", cex = 0.4, pch = 16,
                     add = TRUE)
      graphics::text(x = usr[1] + 0.05 * xwidth,
                     y = usr[3] + 0.08 * yheight,
                     adj = 0,
                     paste("Selected PIs: ", nrow(pis_sub), sep = ""),
                     cex = 1,
                     col = "#fdcdac")
    }
  }

  # plot reference data
  graphics::plot(sf::st_geometry(ned_wh_co_moll),
                 col = NA, border = "#000000", lwd = 1, add = TRUE)
  graphics::plot(sf::st_geometry(ned_wh_st_moll),
                 col = NA, border = "#000000", lwd = 0.75, add = TRUE)

  # label
  usr <- graphics::par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]

  # pds
  if (isTRUE(plot_pds)) {
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.2 * yheight,
                   adj = 0,
                   paste("Available PDs: ", nrow(pds), sep = ""),
                   cex = 1,
                   col = "#1b9377")
    if (!missing(ext)) {
      graphics::text(x = usr[1] + 0.05 * xwidth,
                     y = usr[3] + 0.16 * yheight,
                     adj = 0,
                     paste("Selected PDs: ", nrow(pds_sub), sep = ""),
                     cex = 1,
                     col = "#b3e2cd")
    }
  }

  # pis
  if (isTRUE(plot_pis)) {
    graphics::text(x = usr[1] + 0.05 * xwidth,
                   y = usr[3] + 0.12 * yheight,
                   adj = 0,
                   paste("Available PIs: ", nrow(pis), sep = ""),
                   cex = 1,
                   col = "#d95f02")
    if (!missing(ext)) {
      graphics::text(x = usr[1] + 0.05 * xwidth,
                     y = usr[3] + 0.08 * yheight,
                     adj = 0,
                     paste("Selected PIs: ", nrow(pis_sub), sep = ""),
                     cex = 1,
                     col = "#fdcdac")
    }
  }

  graphics::par(p)
  invisible()
}


#' Calculate and map effective extent of selected centroids
#'
#' The selection of stixel centroids for analysis of PIs and/or PDs yields an
#' effective footprint, or extent, showing the effective location of where the
#' information going into the analysis with PIs and/or PDs is based. While a
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
#' @param pi_pd character; whether to use predictor importance (`"pi"`) or
#'   partial dependence (`"pd"`) for stixel centroids.
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
#' \dontrun{
#'
#' # download and load example data
#' sp_path <- download_data("example_data")
#'
#' # define a spatioremporal extent
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' eff <- calc_effective_extent(path = sp_path, ext = e, pi_pd = "pi")
#'
#' }
calc_effective_extent <- function(path, ext, pi_pd = c("pi", "pd"),
                                  plot = TRUE) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(inherits(ext, "ebirdst_extent"))
  pi_pd <- match.arg(pi_pd)
  stopifnot(is.logical(plot), length(plot) == 1)
  if (all(c(0, 1) == round(ext$t, 2))) {
    warning(paste("Without temporal limits in ext, this function will take",
                  "considerably longer to run and is less informative."))
  }

  # load data
  r_tmplt <- load_raster(product = "template", path = path)
  if (pi_pd == "pi") {
    pipd <- load_pis(path = path)
  } else {
    pipd <- load_pds(path = path)
  }

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
#' @param n integer; the number of colors (≥ 1) to be in the palette.
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
