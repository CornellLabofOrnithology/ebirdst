#' Parse weekly dates from raster layer names
#'
#' The dates corresponding to each layer of a weekly data cube are stored as
#' the layer names. This function converts these to Date objects.
#'
#' @param x [SpatRaster][terra::SpatRaster] object; weekly Status and Trends
#'   data cube.
#'
#' @return Date vector.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # weekly relative abundance
#' abd_weekly <- load_raster(path, "abundance", resolution = "lr")
#'
#' # dates corresponding to each week
#' parse_raster_dates(abd_weekly)
#' }
parse_raster_dates <- function(x) {
  stopifnot(inherits(x, "SpatRaster"))
  dates <- as.Date(names(x), format = "%Y-%m-%d")
  if (any(is.na(dates))) {
    stop("Problem parsing dates stored in input raster layer names. ",
         "Layer names must be formatted 'YYYY-MM-DD'.")
  }
  return(dates)
}

#' Get the Status and Trends week that a date falls into
#'
#' @param dates a vector of dates.
#'
#' @return An integer vector of weeks numbers from 1-52.
#' @export
#' @examples
#' d <- as.Date(c("2016-04-08", "2018-12-31", "2014-01-01", "2018-09-04"))
#' date_to_st_week(d)
date_to_st_week <- function(dates) {
  dv <- seq(from = 0, to = 1, length.out = 52 + 1)
  days <- (as.POSIXlt(dates)$yday + 0.5) / 366

  check_d <- function(x) {
    which(x >= dv[-length(dv)] & x < dv[-1])
  }
  vapply(days, check_d, FUN.VALUE = integer(length = 1))
}


#' Get eBird species code for a set of species
#'
#' Give a vector of species codes, common names, and/or scientific names, return
#' a vector of 6-letter eBird species codes. This function will only look up
#' codes for species for which eBird Status and Trends results exist.
#'
#' @param x character; vector of species codes, common names, and/or scientific
#'   names.
#'
#' @return A character vector of eBird species codes.
#' @export
#'
#' @examples
#' get_species(c("Black-capped Chickadee", "Poecile gambeli", "carchi"))
get_species <- function(x) {
  stopifnot(is.character(x), all(!is.na(x)))
  r <- ebirdst::ebirdst_runs
  x <- tolower(trimws(x))

  # species code
  code <- match(x, tolower(r$species_code))
  # scientific name
  sci <- match(x, tolower(r$scientific_name))
  # common names
  com <- match(x, tolower(r$common_name))
  # combine
  r$species_code[dplyr::coalesce(code, sci, com)]
}


# internal ----

# convert from an iso date to a 0-1 srd date
to_srd_date <- function(x) {
  if (is.character(x)) {
    if (all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x))) {
      x <- as.Date(x)
    } else if (all(grepl("^[0-9]{2}-[0-9]{2}$", x))) {
      # use 2015 since it's not a leap year
      year <- ifelse(x == "02-29", "2016-", "2015-")
      x <- as.Date(paste0(year, x))
    } else {
      stop("Input is not a valid date.")
    }
  } else if (is.numeric(x) || is.integer(x)) {
    x <- as.Date(paste("2016-", x), format = "%Y-%j")
  } else if (!inherits(x, "Date")) {
    stop("Input is not a valid date.")
  }

  # convert to proportion of year
  as.integer(format(x, format = "%j")) / 366
}

# convert from a 0-1 srd date to an iso date
from_srd_date <- function(x, year, iso_date = TRUE) {
  stopifnot(all(x <= 1), all(x >= 0))
  stopifnot(is.logical(iso_date), length(iso_date) == 1)

  if (missing(year)) {
    y <- ebirdst_version()[["version_year"]]
  } else {
    y <- year
  }
  doy <- round(x * 366)
  doy <- pmin(pmax(doy, 1), 365)
  d2015 <- as.Date(x = paste(doy, 2015), format = "%j %Y")
  if (isTRUE(iso_date)) {
    return(as.Date(format(d2015, format = paste0(y, "-%m-%d")),
                   format = "%Y-%m-%d"))
  } else {
    return(format(d2015, format = "%m-%d"))
  }
}

is_integer <- function(x) {
  is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))
}
