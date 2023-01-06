#' Construct a spatiotemporal extent object to subset Status and Trends data
#'
#' `ebirdst_extent` object are used to subset the eBird Status and Trends data
#' spatially and/or temporally. This function constructs these objects.
#'
#' @param x the spatial extent; either a rectangular bounding box (defined as a
#'   vector of numbers representing the coordinates of the boundaries or an
#'   [st_bbox()] object) or a polygon (an [sf] object). See Details for further
#'   explanation of the format of x.
#' @param t the temporal extent; a 2-element vector of the start and end dates
#'   of the temporal extent, provided either as dates (Date objects or strings
#'   in ISO format "YYYY-MM-DD") or numbers between 0 and 1 representing the
#'   fraction of the year. Note that dates can wrap around the year, e.g.
#'   `c("2021-12-01", "2021-01-31") is acceptable. See Details for further
#'   explanation of the format of t. **Leave the argument blank to include the
#'   full year of data.**
#' @param crs coordinate reference system, provided as a `crs` object or
#'   argument to [st_crs()]. Defaults to unprojected, lat/long coordinates (crs
#'   = 4326). **Only required if x is given as a numeric vector defining the
#'   bounding box, ignored in all other cases.**
#' @param ... Additional arguments used by methods.
#'
#' @details The spatial extent, `x`, can be either a rectangular bounding box or
#'   a set of spatial polygons. The bounding box can be defined either as an
#'   [st_bbox()] object or by providing the coordinates of the rectangle edges
#'   directly as a named vector with elements xmin, xmax, ymin, and ymax (note
#'   that latitude and longitude correspond to y and x, respectively). In this
#'   latter case, a coordinate reference system must be provided explicitly via
#'   the `crs` argument (`crs = 4326` is the default and is a short form for
#'   unprojected lat/long coordinates). For a polygon spatial extent, `x` should
#'   be either an [sf] or [sfc] object (with feature type `POLYGON` or
#'   `MULTIPOLYGON`) from the `sf` package. To import data from a Shapefile or
#'   GeoPackage into this format, use [read_sf()].
#'
#'   The temporal extent defines the start and end dates of the time period.
#'   These are most easily provided as Date objects or date strings in ISO
#'   format ("YYYY-MM-DD"). If dates are defined as strings, the year can be
#'   omitted (i.e. "MM-DD"). Alternatively, dates can be defined in terms of
#'   fractions of the year, e.g. `t = c(0.25, 0.5) ` would subset to data within
#'   the second quarter of the year. In all cases, dates can wrap around the
#'   year, e.g. c("2021-12-01", "2021-01-31") would subset to data in December
#'   or January.
#'
#' @return An `ebirdst_extent` object consisting of a list with three elements:
#'   the spatial extent `extent`, the temporal extent `t`, and `type` (either
#'   "bbox" or "polygon").
#'
#' @export
#' @examples
#' # bounding box of the north eastern united stats as a numeric vector
#' bb_vec <- c(xmin = -80, xmax = -70, ymin = 40, ymax = 47)
#' ebirdst_extent(bb_vec)
#'
#' # bbox object
#' bb <- sf::st_bbox(bb_vec, crs = 4326)
#' ebirdst_extent(bb)
#'
#' # polygon imported from a shapefile
#' poly <- sf::read_sf(system.file("shape/nc.shp", package="sf"))
#' ebirdst_extent(poly)
#'
#' # subset to january
#' ebirdst_extent(bb, t = c("2021-01-01", "2021-01-31"))
#'
#' # dates can wrap around, e.g. to use dec-jan
#' ebirdst_extent(bb, t = c("2021-12-01", "2021-01-31"))
#'
#' # dates can also be given without an associated year
#' ebirdst_extent(bb, t = c("12-01", "01-31"))
ebirdst_extent <- function(x, t, ...) {
  UseMethod("ebirdst_extent")
}

#' @export
#' @describeIn ebirdst_extent bounding box created with [st_bbox()]
ebirdst_extent.bbox <- function(x, t, ...) {
  if (is.na(sf::st_crs(x)$proj4string)) {
    stop("A CRS must be defined for the spatial extent x.")
  }
  # default to full year
  if (missing(t)) {
    t <- c(0, 1)
  }

  # convert temporal extent
  t <- process_t_extent(t)

  structure(
    list(type = "bbox", extent = x, t = t),
    class = "ebirdst_extent"
  )
}

#' @export
#' @describeIn ebirdst_extent bounding box given as edges
ebirdst_extent.numeric <- function(x, t, crs = 4326, ...) {
  stopifnot(is.numeric(x), length(x) == 4, all(!is.na(x)),
            !is.null(names(x)),
            all(c("xmin", "ymin", "xmax", "ymax") %in% names(x)),
            x["xmin"] < x["xmax"], x["ymin"] < x["ymax"])

  ebirdst_extent.bbox(x = sf::st_bbox(obj = x, crs = crs), t = t)
}

#' @export
#' @describeIn ebirdst_extent polygons as [sfc] spatial feature column
ebirdst_extent.sfc <- function(x, t, ...) {
  if (any(!sf::st_is(x, c("MULTIPOLYGON", "POLYGON")))) {
    stop("Spatial extent must consist of polygon features.")
  }
  if (all(is.na(sf::st_crs(x)$proj4string))) {
    stop("A CRS must be defined for the spatial extent x.")
  }
  # default to full year
  if (missing(t)) {
    t <- c(0, 1)
  }

  # convert temporal extent
  t <- process_t_extent(t)

  structure(
    list(type = "polygon", extent = sf::st_combine(x), t = t),
    class = "ebirdst_extent"
  )
}


#' @export
#' @describeIn ebirdst_extent polygons as [sf] object
ebirdst_extent.sf <- function(x, t, ...) {
  ebirdst_extent.sfc(x = sf::st_geometry(x), t = t)
}


#' Transform a spatiotemporal extent to a different CRS
#'
#' Transform an eBird Status and Trends extent object to a different
#' coordinate reference system. This is most commonly required to transform the
#' extent to the sinusoidal CRS used by the eBird Status and Trends rasters.
#'
#' @param x [ebirdst_extent] object; a spatiotemporal extent.
#' @param crs coordinate references system, given either as a proj4string, an
#'   integer EPSG code, or a `crs` object generated with [st_crs()].
#'
#' @return An [ebirdst_extent] object in the new CRS.
#' @export
#'
#' @examples
#' # construct an ebirdst_extent object
#' bb_vec <- c(xmin = -80, xmax = -70, ymin = 40, ymax = 47)
#' bb <- sf::st_bbox(bb_vec, crs = 4326)
#' bb_ext <- ebirdst_extent(bb)
#'
#' # transform to sinusoidal projection of rasters
#' sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
#' project_extent(bb_ext, crs = sinu)
#'
#' # also works on polygon extents
#' poly <- sf::read_sf(system.file("shape/nc.shp", package="sf"))
#' poly_ext <- ebirdst_extent(poly)
#' project_extent(poly_ext, crs = sinu)
project_extent <- function(x, crs) {
  stopifnot(inherits(x, "ebirdst_extent"))

  if (x$type == "bbox") {
    bb_poly <- sf::st_as_sfc(x$extent)
    bb_poly_proj <- sf::st_transform(bb_poly, crs = crs)
    x$extent <- sf::st_bbox(bb_poly_proj)
  } else if (x$type == "polygon") {
    x$extent <- sf::st_transform(x$extent, crs = crs)
  } else {
    stop("Invalid ebirst_extent object.")
  }
  return(x)
}


#' @export
print.ebirdst_extent <- function(x, ...) {
  cat("eBird Status & Trends extent: \n")

  # spatial
  if (x$type == "polygon") {
    cls <- substr(class(x$extent)[1], 5, nchar(class(x$extent)[1]))
    cat(paste("  Polygon:", length(x$extent), cls, "features\n"))
  }
  bb <- signif(sf::st_bbox(x$extent), options("digits")$digits)
  cat("  Bounding box: ")
  cat(paste(paste(names(bb), bb[]), collapse = "; "))

  cat("\n")
  cat(paste("  CRS:", sf::st_crs(x$extent)$proj4string, "\n"))

  # temporal
  t_date <- format(from_srd_date(x$t), format = "%m-%d")
  cat(paste("  Temporal extent:", t_date[1], " - ", t_date[2], "\n"))
  invisible(x)
}


# internal ----
process_t_extent <- function(t) {
  stopifnot(length(t) == 2, all(!is.na(t)))

  if (is.numeric(t)) {
    stopifnot(all(t <= 1), all(t >= 0))
  } else if (inherits(t, "Date") || is.character(t)) {
    t <- to_srd_date(t)
  } else {
    stop("Unrecognized format for temporal extent t.")
  }
  return(t)
}
