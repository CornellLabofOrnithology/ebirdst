#' Subset eBird Status and Trends data spatiotemporally
#'
#' Spatiotemporally subset the raster or tabular eBird Status and Trends data.
#' The spatiotemporal extent should be defined using [ebirdst_extent()].
#'
#' @param x eBird Status and Trends data to subset; either a Raster object with
#'   52 layers (one for each week) or a `data.frame` with PI or PD data.
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to
#'   filter the data to.
#'
#' @return eBird Status and Trends data in the same format as the input data.
#' @export
#'
#' @examples \dontrun{
#' # bbox for southern michigan in may
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' # download and load example data
#' sp_path <- ebirdst_download("example_data")
#' pis <- load_pis(sp_path)
#' abd <- load_raster(product = "abundance", sp_path)
#'
#' # subset
#' abd_ss <- ebirdst_subset(abd, ext = e)
#' pis_ss <- ebirdst_subset(pis, ext = e)
#' }
ebirdst_subset <- function(x, ext) {
  UseMethod("ebirdst_subset")
}

#' @export
#' @describeIn ebirdst_subset PI or PD data
ebirdst_subset.data.frame <- function(x, ext) {
  stopifnot(all(c("date", "lon", "lat") %in% names(x)))
  stopifnot(all(x$date >= 0), all(x$date <= 1))
  stopifnot(inherits(ext, "ebirdst_extent"))

  # temporal filtering
  if (!identical(ext$t, c(0, 1))) {
    if (ext$t[1] <= ext$t[2]) {
      x <- x[x$date > ext$t[1] & x$date <= ext$t[2], ]
    } else {
      x <- x[x$date > ext$t[1] | x$date <= ext$t[2], ]
    }
  } else {
    warning("Temporal extent missing or incomplete. This may be intentional.")
  }

  # spatial filtering
  e_ll <- project_extent(ext, crs = 4326)
  if (ext$type == "bbox") {
    b <- e_ll$extent
    x <- x[x$lon > b["xmin"] & x$lon <= b["xmax"] &
             x$lat > b["ymin"] & x$lat <= b["ymax"], ]
  } else if (ext$type == "polygon") {
    x_sf <- sf::st_as_sf(x, coords = c("lon", "lat"), crs = 4326)
    is_in <- suppressMessages(
      sf::st_intersects(x_sf, e_ll$extent, sparse = FALSE)
    )
    if (!is.matrix(is_in) || ncol(is_in) != 1) {
      stop("Problem with ebirdst_extent object.")
    }
    x <- x[is_in[, 1, drop = TRUE], ]
  } else {
    stop("Invalid ebirdst_extent object.")
  }
  return(x)
}


#' @export
#' @describeIn ebirdst_subset  PI or PD data as an `sf` object
ebirdst_subset.sf <- function(x, ext) {
  stopifnot("date" %in% names(x))
  stopifnot(all(x$date >= 0), all(x$date <= 1))
  stopifnot(inherits(ext, "ebirdst_extent"))

  # temporal filtering
  if (!identical(ext$t, c(0, 1))) {
    if (ext$t[1] <= ext$t[2]) {
      x <- x[x$date > ext$t[1] & x$date <= ext$t[2], ]
    } else {
      x <- x[x$date > ext$t[1] | x$date <= ext$t[2], ]
    }
  } else {
    warning("Temporal extent missing or incomplete. This may be intentional.")
  }

  # spatial filtering
  e_ll <- project_extent(ext, crs = sf::st_crs(x))
  if (ext$type == "bbox") {
    b <- sf::st_as_sfc(e_ll$extent)
  } else if (ext$type == "polygon") {
    b <- e_ll$extent
  } else {
    stop("Invalid ebirdst_extent object.")
  }
  is_in <- suppressMessages(
    sf::st_intersects(x, b, sparse = FALSE)
  )
  if (!is.matrix(is_in) || ncol(is_in) != 1) {
    stop("Problem with ebirdst_extent object.")
  }
  x <- x[is_in[, 1, drop = TRUE], ]
}

#' @export
#' @describeIn ebirdst_subset Status and Trends rasters
ebirdst_subset.Raster <- function(x, ext) {
  if((raster::nlayers(x) != 52)) {
    stop(paste("The input Raster object must be full stack or brick of 52",
               "layers as originally provided."))
  }
  stopifnot(inherits(ext, "ebirdst_extent"))

  # temporal filtering
  x <- label_raster_stack(x)
  if (!identical(ext$t, c(0, 1))) {
    r_dates <- to_srd_date(parse_raster_dates(x))
    if (ext$t[1] <= ext$t[2]) {
      r_in <- r_dates > ext$t[1] & r_dates <= ext$t[2]
    } else {
      r_in <- r_dates > ext$t[1] | r_dates <= ext$t[2]
    }
    x <- x[[which(r_in)]]
  } else {
    warning("Temporal extent missing or incomplete. This may be intentional.")
  }

  # spatial filtering
  e_ll <- project_extent(ext, crs = sf::st_crs(x))
  if (ext$type == "bbox") {
    x <- raster::crop(x, bbox_to_extent(e_ll$extent))
  } else if (ext$type == "polygon") {
    x <- raster::trim(
      raster::mask(
        raster::crop(x, bbox_to_extent(sf::st_bbox(e_ll$extent))),
        sf::st_sf(e_ll$extent)),
      values = NA)
  } else {
    stop("Invalid ebirdst_extent object.")
  }
  return(x)
}


bbox_to_extent <- function(x) {
  raster::extent(x[["xmin"]], x[["xmax"]], x[["ymin"]], x[["ymax"]])
}
