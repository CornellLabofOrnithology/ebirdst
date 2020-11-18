


#' Label data cubes with the week date for each band
#'
#' The data cubes are saved as GeoTIFFs, which don't allow for band labels. For
#' convenience, this function labels the layers of a data cube once it has been
#' loaded with the week dates for each band.
#'
#' @param x `RasterStack` or `RasterBrick`; original eBird Status and Trends
#'   data cube with 52 bands, one for each week.
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
label_raster_stack <- function(x) {
  stopifnot(inherits(x, "Raster"))

  # check length
  if((raster::nlayers(x) != 52)) {
    stop(paste("The input Raster* object must be a full cube of 52",
               "layers as originally provided."))
  }

  srd_date_vec <- seq(from = 0, to = 1, length = 52 + 1)
  srd_date_vec <- (srd_date_vec[1:52] + srd_date_vec[2:(52 + 1)]) / 2
  srd_date_vec <- round(srd_date_vec, digits = 4)

  year_seq <- 2015
  p_time <- strptime(x = paste(round(srd_date_vec * 366), year_seq), "%j %Y")
  date_names <- paste(paste0("w", ebirdst_version()["data_version"]),
                      formatC(p_time$mon + 1, width = 2, format = "d",
                              flag = "0"),
                      formatC(p_time$mday, width = 2, format = "d",
                              flag = "0"),
                      sep = ".")

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
#' # download and load example abundance data
#' sp_path <- ebirdst_download("example_data")
#' abd <- load_raster(sp_path, "abundance")
#'
#' # parse dates
#' parse_raster_dates(abd)
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
  dv <- seq(from = 0, to = 1, length = 52 + 1)
  days <- (as.POSIXlt(dates)$yday + 0.5) / 366

  check_d <- function(x) {
    which(x >= dv[-length(dv)] & x <= dv[-1])
  }
  vapply(days, check_d, FUN.VALUE = integer(length = 1))
}


ebirdst_version <- function() {
  c(data_version = 2019, release_year = 2020)
}
