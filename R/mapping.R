#' Calculate the spatial extent of non-zero data in a raster
#'
#' eBird Status and Trends data cubes are defined over broad areas, filling in
#' regions where the species doesn't occur with zeros (predicted absences) or
#' NAs (regions where models weren't fit). When producing maps, it's best to
#' only display the spatial extent where the species occurs. To show determine
#' an ideal extent for mapping, this function trims away 0 and NA values. When
#' called on a `RasterStack` (e.g., a data cube consisting of all 52 weeks),
#' this function returns the extent of occurrence across all layers.To access a
#' pre-calculated extent for the full annual cycle use
#' [load_fac_map_parameters()].
#'
#' @param x `Raster` object; either a full 52-week data cube or a subset.
#' @param aggregate logical; whether data should be aggregated by a factor of 3
#'   in each dimension prior to calculating the extent. When working with the
#'   high resolution cubes, data should be aggregated otherwise processing times
#'   will be extremely long.
#'
#' @return The extent of occurrence as a `raster` [Extent] object.
#'
#' @export
#'
#' @examples
#' # simple toy example
#' r <- raster::raster(nrow = 100, ncol = 100)
#' r[5025:5075] <- 1
#' raster::extent(r)
#' calc_full_extent(r)
#'
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load abundance data
#' abd <- load_raster(path, "abundance")
#'
#' # calculate full extent
#' map_extent <- calc_full_extent(abd)
#'
#' # plot
#' raster::plot(abd[[20]], axes = FALSE, ext = map_extent)
#' }
calc_full_extent <- function(x, aggregate = TRUE) {
  stopifnot(inherits(x, "Raster"))
  stopifnot(is.logical(aggregate), length(aggregate) == 1)

  # aggregate stack for speed, otherwise everything else takes too long
  e_input <- raster::extent(x)
  if (isTRUE(aggregate)) {
    x <- raster::aggregate(x, fact = 3)
  }

  # convert 0s to NAs, otherwise trimming is slow and the extent is too broad
  x[x == 0] <- NA

  # work on individual layers because trim has a bug when applied to stacks
  # also this approach is faster
  e <- c(NA_real_, NA_real_, NA_real_, NA_real_)
  for (i in seq_len(raster::nlayers(x))) {
    e_trim <- raster::extent(raster::trim(x[[i]]))
    e[1] <- min(e[1], e_trim[1], na.rm = TRUE)
    e[2] <- max(e[2], e_trim[2], na.rm = TRUE)
    e[3] <- min(e[3], e_trim[3], na.rm = TRUE)
    e[4] <- max(e[4], e_trim[4], na.rm = TRUE)
  }
  e <- raster::extent(e)

  # sometimes extent calculations get weird and you'll get a very broad
  # extent that goes further than you want, so check against the input

  e[1] <- max(e[1], e_input[1], na.rm = TRUE)
  e[2] <- min(e[2], e_input[2], na.rm = TRUE)
  e[3] <- max(e[3], e_input[3], na.rm = TRUE)
  e[4] <- min(e[4], e_input[4], na.rm = TRUE)

  return(e)
}


#' Calculates relative abundance bins (breaks) based for mapping
#'
#' Mapping species abundance across the full-annual cycle presents a challenge,
#' in that patterns of concentration and dispersion in abundance change
#' throughout the year, making it difficult to define color bins that suit all
#' seasons and accurately reflect the detail of abundance predictions. To
#' address this, when mapping the relative abundance data, we recommend using
#' quantile bins based on the underlying count distribution, adjusted
#' according to the relative abundance distribution. To access pre-calculated
#' bins for the full annual cycle use [load_fac_map_parameters()].
#'
#' @param abundance `Raster` object; eBird Status and Trends relative abundance
#'   data cube, or a subset of the data cube.
#' @param count `Raster` object; eBird Status and Trends count data cube, or a
#'   subset of the data cube.
#'
#' @return A numeric vector defining the breaks of the relative abundance bins.
#'   In addition, the `labels` attribute of this vector provides the 5th, 50th,
#'   and 95th quantile of abundance, which can be used to label the bottom,
#'   middle, and top of a legend.
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
#' # abundance data
#' abd <- load_raster(path, "abundance")
#' # count data
#' cnt <- load_raster(path, "count")
#'
#' # calculate bins for a single week for this example
#' bins <- calc_bins(abd, cnt)
#' }
calc_bins <- function(abundance, count) {
  stopifnot(inherits(abundance, "Raster"))
  stopifnot(inherits(count, "Raster"))

  if (all(is.na(suppressWarnings(raster::maxValue(abundance)))) &&
      all(is.na(suppressWarnings(raster::minValue(abundance)))) &&
      all(is.na(suppressWarnings(raster::maxValue(count)))) &&
      all(is.na(suppressWarnings(raster::minValue(count))))) {
    stop("Input Raster* objects must have non-NA values.")
  }
  if (all(raster::maxValue(abundance) == 0) &&
      all(raster::maxValue(count) == 0)) {
    stop("Raster must have at least 2 non-zero values to calculate bins")
  }

  # abundance
  v <- as.vector(raster::getValues(abundance))
  v <- as.numeric(stats::na.omit(v))
  v <- v[v > 0]
  if (length(v) <= 1) {
    stop("Raster must have at least 2 non-zero values to calculate bins")
  }
  abd_rng <- range(v, na.rm = TRUE)
  abd_5th <- stats::quantile(v, 0.05)
  rm(v)

  # count
  v <- as.vector(raster::getValues(count))
  v <- as.numeric(stats::na.omit(v))
  v <- v[v > 0]
  # quantile bins based on counts
  b <- stats::quantile(v, probs = seq(0, 1, by = 0.05), na.rm = TRUE)

  # adjust for abundance
  b[1] <- abd_rng[1]
  b <- b[b <= abd_rng[2]]
  b <- c(b, abd_rng[2])
  if (b[2] > abd_5th) {
    b <- sort(c(b, abd_5th))
  }
  b <- unname(b)
  # labels
  attr(b, "labels") <- c(b[2], stats::median(b), b[length(b) - 1])

  return(b)
}


#' eBird Status and Trends color palettes for mapping
#'
#' Generate the color palettes used for the eBird Status and Trends relative
#' abundance maps.
#'
#' @param n integer; the number of colors to be in the palette.
#' @param season character; the season to generate colors for or "weekly" to
#'   get the color palette used in the weekly abundance animations.
#'
#' @return A character vector of hex color codes.
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
  col_zero <- "#e6e6e6"
  if (season == "weekly") {
    plsm <- rev(viridisLite::plasma(n - 1, end = 0.9))
    plsm <- stringr::str_remove(plsm, "FF$")
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
