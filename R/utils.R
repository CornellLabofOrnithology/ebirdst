#' Label data cubes with the week date for each band
#'
#' The data cubes are saved as GeoTIFFs, which don't allow for band labels. For
#' convenience, this function labels the layers of a data cube once it has been
#' loaded with the week dates for each band.
#'
#' @param x `RasterStack` or `RasterBrick`; eBird Status and Trends data cube,
#'   typically with 52 bands, one for each week.
#' @param weeks vector of dates corresponding to the weeks of each layer in `x`.
#'   This argument should be used only rarely, when you're labelling a cube with
#'   fewer that the usual 52 weeks of predictions.
#'
#' @return A `RasterStack` or `RasterBrick` with names assigned for the dates in
#'   the format of "wYYYY.MM.DD" per raster package constraints. The `Raster*`
#'   objects do not allow the names to start with a number, nor are they allowed
#'   to contain "-", so it is not possible to store the date in an ISO compliant
#'   format. Use `parse_raster_dates()` to convert the layer names to dates.
#'
#' @export
#'
#' @examples
#' # download and load example abundance data
#' sp_path <- ebirdst_download("example_data")
#' abd <- load_raster(sp_path, "abundance")
#'
#' # label
#' abd <- label_raster_stack(abd)
#' names(abd)
label_raster_stack <- function(x, weeks = NULL) {
  stopifnot(inherits(x, "Raster"))

  if (is.null(weeks)) {
    if((raster::nlayers(x) != 52)) {
      stop("The input Raster* object must be a full cube of 52 weeks unless ",
           "a weeks argument is provided to label_raster_stack().")
    }
    weeks <- ebirdst::ebirdst_weeks$date
  }
  stopifnot(inherits(weeks, "Date"))

  # check lengths
  if (raster::nlayers(x) != length(weeks)) {
    stop("The number of raster layers must match the number of weeks.")
  }

  date_names <- format(weeks, "w%Y.%m.%d")
  names(x) <- date_names

  return(x)
}


#' Parse data cube layer names into dates
#'
#' [label_raster_stack()] labels the layers of a data cube with the associated
#' week dates in the format of "wYYYY.MM.DD", because of constraints in the
#' `raster` package. This function converts that character vector into an ISO
#' compliant Date vector.
#'
#' @param x `Raster` object; labeled Status and Trends data cube.
#'
#' @return Date vector.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download and load example abundance data
#' sp_path <- ebirdst_download("example_data")
#' abd <- load_raster(sp_path, "abundance")
#'
#' # parse dates
#' parse_raster_dates(abd)
#' }
parse_raster_dates <- function(x) {
  stopifnot(inherits(x, "Raster"))
  if (!all(grepl("w[0-9]{4}\\.[0-9]{2}\\.[0-9]{2}", names(x)))) {
    stop("Raster names not in correct format, call label_raster_stack() first.")
  }

  as.Date(names(x), format = "w%Y.%m.%d")
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


# internal functions ----

ebirdst_version <- function() {
  c(data_version = 2020, release_year = 2021)
}

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
    y <- ebirdst_version()["data_version"]
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
