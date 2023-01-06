#' Load eBird Status and Trends raster data cubes
#'
#' Each of the eBird Status and Trends raster products is packaged as a GeoTIFF
#' file representing predictions on a regular grid. The core products are
#' occurrence, count, relative abundance, and percent of population. This
#' function loads one of the available data products into R as a
#' [SpatRaster][terra::SpatRaster] object.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#' @param product character; Status and Trends product to load: occurrence,
#'   count, relative abundance, or percent of population. See Details for a
#'   detailed explanation of each of these products.
#' @param period character; temporal period of the estimation. The Status and
#'   Trends models make predictions for each week of the year; however, as a
#'   convenience, data are also provided summarized at the seasonal or annual
#'   ("full-year") level.
#' @param metric character; by default, the weekly products provide estimates of
#'   the median value (`metric = "median"`) and the summarized products are the
#'   cell-wise mean across the weeks within the season (`metric = "mean"`).
#'   However, additional variants exist for some of the products. For the weekly
#'   relative abundance, confidence intervals are provided: specify `metric =
#'   "lower"` to get the 10th quantile or `metric = "upper"` to get the 90th
#'   quantile. For the seasonal and annual products, the cell-wise maximum
#'   values across weeks can be obtained with `metric = "max"`.
#' @param resolution character; the resolution of the raster data to load. The
#'   default is to load the native ~3 km resolution (`"hr"`); however, for some
#'   applications 9 km (`"mr"`) or 27 km (`"lr"`) data may be suitable.
#'
#' @details The core Status and Trends data products provide weekly estimates
#'   across a regular spatial grid. They are packaged as rasters with 52 layers,
#'   each corresponding to estimates for a week of the year, and we refer to
#'   them as "cubes" (e.g. the "relative abundance cube"). All estimates are the
#'   median expected value for a standard 1km, 1 hour eBird Traveling Count by
#'   an expert eBird observer at the optimal time of day and for optimal weather
#'   conditions to observe the given species. These products are:
#'
#' - `occurrence`: the expected probability (0-1) of occurrence a species.
#' - `count`: the expected count of a species, conditional on its occurrence at
#' the given location.
#' - `abundance`: the expected relative abundance of a species, computed as the product of
#' the probability of occurrence and the count conditional on occurrence.
#' - `percent-population`: the percent of the total relative abundance within
#' each cell. This is a derived product calculated by dividing each cell value
#' in the relative abundance raster with the total abundance summed across all
#' cells.
#'
#' In addition to these weekly data cubes, this function provides access to data
#' summarized over different periods. Seasonal cubes are produced by taking the
#' cell-wise mean or max across the weeks within each season. The boundary dates
#' for each season are species specific and are available in `ebirdst_runs`, and
#' if a season failed review no associated layer will be included in the cube.
#' In addition, full-year summaries provide the mean or max across all weeks of
#' the year that fall within a season that passed review. Note that this is not
#' necessarily all 52 weeks of the year. For example, if the estimates for the
#' non-breeding season failed expert review for a given species, the full-year
#' summary for that species will not include the weeks that would fall within
#' the non-breeding season.
#'
#' @return For the weekly cubes, a [SpatRaster][terra::SpatRaster] with 52
#'   layers for the given product, labeled by week. Seasonal cubes will have up
#'   to four layers labeled according to the seasons. The full-year products
#'   will have a single layer.
#'
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
#' # note that only low resolution (lr) data are available for the example data
#' abd_weekly <- load_raster(path, "abundance", resolution = "lr")
#'
#' # the weeks for each layer are stored in the layer names
#' names(abd_weekly)
#' # and can be converted to Date objects with
#' parse_raster_dates(abd_weekly)
#'
#' # max seasonal abundance
#' abd_seasonal <- load_raster(path, "abundance",
#'                             period = "seasonal", metric = "max",
#'                             resolution = "lr")
#' # available seasons in stack
#' names(abd_seasonal)
#' # subset to just breeding season abundance
#' abd_seasonal[["breeding"]]
#' }
load_raster <- function(path,
                        product = c("abundance",
                                    "count",
                                    "occurrence",
                                    "percent-population"),
                        period = c("weekly",
                                   "seasonal",
                                   "full-year"),
                        metric = NULL,
                        resolution = c("hr", "mr", "lr")) {

  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  product <- match.arg(product)
  period <- match.arg(period)
  resolution <- match.arg(resolution)

  # check that the geotiff driver is installed
  drv <- terra::gdal(drivers = TRUE)
  drv <- drv$name[stringr::str_detect(drv$can, "read")]
  if (!"GTiff" %in% drv) {
    stop("GDAL does not have GeoTIFF support. GeoTIFF support is required to ",
         "load Status and Trends raster data.")
  }

  # load config file
  p <- load_config(path)
  species_code <- p[["species_code"]]
  v <- ebirdst_version()[["version_year"]]
  is_example <- stringr::str_detect(species_code, "-example")

  if (is_example && !resolution == "lr") {
    stop("The example data only contains low-resolution (lr) estimates.")
  }

  # full year products only available for migrants
  if (p$is_resident && period == "full-year") {
    stop("Full-year products are not available for residents, use ",
         "period = 'seasonal' instead.")
  }

  # construct file name and path
  if (period == "weekly") {
    # assess which metric is being requested
    if (is.null(metric)) {
      metric <- "median"
    }
    if (product == "abundance") {
      if (!metric %in% c("median", "lower", "upper")) {
        stop("Valid metrics for weekly abundance data are 'median', 'lower', ",
             "or 'upper'")
      }
    } else {
      if (metric != "median") {
        stop("For this product, metric must be 'median'")
      }
    }

    # construct filename
    file <- stringr::str_glue("{species_code}_{product}_{metric}",
                              "_{resolution}_{v}.tif")
    file <- file.path(path, "weekly", file)
  } else {
    # assess which metric is being requested
    if (is.null(metric)) {
      metric <- "mean"
    }
    if (!metric %in% c("mean", "max")) {
      stop("Valid metrics for seasonal or full-year data are 'mean' or 'max.'")
    }

    # construct filename
    file <- stringr::str_glue("{species_code}_{product}_{period}_{metric}",
                              "_{resolution}_{v}.tif")
    file <- file.path(path, "seasonal", file)
  }

  # check existence of target file
  if (!file.exists(file)) {
    stop("The file for the requested product does not exist: \n  ", file)
  }

  # load and return raster stack
  return(terra::rast(file))
}


#' Load seasonal eBird Status and Trends range polygons
#'
#' Range polygons are defined as the boundaries of non-zero seasonal relative
#' abundance estimates, which are then (optionally) smoothed to produce more
#' aesthetically pleasing polygons using the `smoothr` package.
#'
#' @inheritParams load_raster
#' @param resolution character; the raster resolution from which the range
#'   polygons were derived.
#' @param smoothed logical; whether smoothed or unsmoothed ranges should be
#'   loaded.
#'
#' @return An `sf` update containing the seasonal range boundaries, with each
#'   season provided as a different feature.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load smoothed ranges
#' # note that only low res data are provided for the example data
#' ranges <- load_range(path, resolution = "lr")
#' }
load_ranges <- function(path, resolution = c("mr", "lr"), smoothed = TRUE) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(is.logical(smoothed), length(smoothed) == 1)
  resolution <- match.arg(resolution)

  # load config file
  p <- load_config(path)
  species_code <- p[["species_code"]]
  v <- ebirdst_version()[["version_year"]]
  is_example <- stringr::str_detect(species_code, "-example")

  if (is_example && !resolution == "lr") {
    stop("The example data only contains low-resolution (lr) estimates.")
  }

  # define filename
  label <- ifelse(smoothed, "smooth", "raw")
  file <- stringr::str_glue("{species_code}_range_{label}",
                            "_{resolution}_{v}.gpkg")
  file <- file.path(path, "ranges", file)

  # check existence of target file
  if (!file.exists(file)) {
    stop("The file for the requested product does not exist: \n  ", file)
  }

  # load polygons
  p <- sf::read_sf(dsn = file, layer = "range")

  return(p)
}


#' Load eBird Status and Trends predictor importance data
#'
#' Loads the predictor importance (PI) data from the stixel_summary.db sqlite
#' database. PI estimates are provided for each stixel over which a model was
#' run and are identified by a unique stixel ID in addition to the coordinates
#' of the stixel centroid.
#'
#' @inheritParams load_raster
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to filter the
#'   data to. The spatial component of the extent object must be provided in
#'   unprojected, latitude-longitude coordinates.
#' @param model character; whether to make estimates for the occurrence or count
#'   model.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return Data frame, or [sf] object if `return_sf = TRUE`, containing PI
#'   estimates for each stixel for either the occurrence or count models. The
#'   data are provided in a 'wide' format, with each row corresponding to the PI
#'   estimates for a give stixel for the occurrence count model, and the
#'   relative importance of each predictor in columns. Stixels are identified by
#'   a unique `stixel_id`, and the centroid of the stixel in space and time is
#'   specified by the `latitude`, `longitude`, and `day_of_year` columns. The
#'   column `predictor` provides a code specifying the predictor variable. These
#'   codes can be looked up in `ebirdst_predictors` for a brief description.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load predictor importance for the occurrence model
#' pis <- load_pis(path)
#'
#' # plot the top 15 predictor importances
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pis(pis, ext = e, n_top_pred = 15, by_cover_class = TRUE)
#' }
load_pis <- function(path, ext, model = c("occurrence", "count"),
                     return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }
  model <- match.arg(model)
  table <- paste0(model, "_pis")

  db_file <- file.path(path, "stixel_summary.db")
  if(!file.exists(db_file)) {
    stop("The file 'stixel_summary.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  if (missing(ext)) {
    sql <- stringr::str_glue("SELECT p.*, s.latitude, s.longitude, s.day_of_year ",
                             "FROM {table} AS p INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id;")
  } else {
    sql <- stringr::str_glue("SELECT p.*, s.latitude, s.longitude, s.day_of_year ",
                             "FROM {table} AS p ",
                             "INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id ",
                             "{sql_extent_subset(ext)};")
  }

  # extract from database
  pis <- DBI::dbGetQuery(db, sql)
  pis <- dplyr::tibble(pis)
  DBI::dbDisconnect(db)

  # clean up names
  pis <- pis[, c("stixel_id", "latitude", "longitude", "day_of_year",
                 "covariate", "pi")]
  pis <- dplyr::rename(pis, predictor = "covariate", importance = "pi")
  pis[["predictor"]] <- transform_predictor_names(pis[["predictor"]])

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(pis[, c("latitude", "longitude",
                                                "day_of_year")])
  if (any(!has_centroid)) {
    warning("Removing ", sum(!has_centroid),
            " stixels with missing centroids data.")
    pis <- pis[has_centroid, ]
  }

  # subset to polygon if one was provided
  if (!missing(ext) && ext$type == "polygon") {
    pis <- ebirdst_subset(pis, ext = ext)
  }

  if (isTRUE(return_sf)) {
    pis <- sf::st_as_sf(pis, coords = c("longitude", "latitude"), crs = 4326)
  }

  attr(pis, "model") <- model
  return(pis)
}


#' Load eBird Status and Trends partial dependence data
#'
#' Partial dependence (PD) plots depict the relationship between the modeled
#' occurrence probability and each of the predictor variables used in the model.
#' Status and Trends provides the data to generate these plots for every stixel.
#'
#' @inheritParams load_pis
#'
#' @return Data frame, or [sf] object if `return_sf = TRUE`, containing PD
#'   estimates for each stixel for either the occurrence or count model. The
#'   data frame will have the following columns:
#'   - `stixel_id`: unique stixel identifier
#'   - `latitude` and `longitude`: stixel centroid
#'   - `day_of_year`: center day of year for stixel
#'   - `predictor`: name of the predictor that the PD data correspond to, for a
#'   full list of predictors consult the [ebirdst_predictors] data frame
#'   - `predictor_value`: value of the predictor variable at which PD is
#'   evaluated
#'   - `response`: predicted response, occurrence or count, at the given value
#'   of the predictor averaged across all the values of the other predictors
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load partial dependence data for occurrence model
#' pds <- load_pds(path)
#'
#' # plot the top 15 predictor importances
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pds(pds, "solar_noon_diff_mid", ext = e, n_bs = 5)
#' }
load_pds <- function(path, ext, model = c("occurrence", "count"),
                     return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }
  model <- match.arg(model)
  table <- paste0(model, "_pds")

  db_file <- file.path(path, "stixel_summary.db")
  if(!file.exists(db_file)) {
    stop("The file 'stixel_summary.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  if (missing(ext)) {
    sql <- stringr::str_glue("SELECT p.*, s.latitude, s.longitude, s.day_of_year ",
                             "FROM {table} AS p INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id;")
  } else {
    sql <- stringr::str_glue("SELECT p.*, s.latitude, s.longitude, s.day_of_year ",
                             "FROM {table} AS p ",
                             "INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id ",
                             "{sql_extent_subset(ext)};")
  }

  # extract from database
  pds <- DBI::dbGetQuery(db, sql)
  pds <- dplyr::tibble(pds)
  DBI::dbDisconnect(db)

  # clean up names
  pds <- pds[, c("stixel_id", "latitude", "longitude", "day_of_year",
                 "covariate", "predictor_value", "response")]
  pds <- dplyr::rename(pds, predictor = "covariate")
  pds[["predictor"]] <- transform_predictor_names(pds[["predictor"]])

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(pds[, c("latitude", "longitude",
                                                "day_of_year")])
  if (any(!has_centroid)) {
    warning("Removing ", sum(!has_centroid),
            " stixels with missing centroids data.")
    pds <- pds[has_centroid, ]
  }

  # subset to polygon if one was provided
  if (!missing(ext) && ext$type == "polygon") {
    pds <- ebirdst_subset(pds, ext = ext)
  }

  if (isTRUE(return_sf)) {
    pds <- sf::st_as_sf(pds, coords = c("longitude", "latitude"), crs = 4326)
  }

  attr(pds, "model") <- model
  return(pds)
}


#' Load summary data for eBird Status and Trends stixels
#'
#' eBird Status and Trends divides space and time into variably sized "stixels"
#' within which individual base models are fit. The process of stixelization is
#' performed many times and the prediction at any given point is the median of
#' the predictions from all the stixels that that point falls in. This function
#' loads summary statistics for each stixel, for example, the size of the
#' stixels, the number of observations within each stixel, and a suite of
#' predictive performance metrics (PPMs) for the model fit within that stixel.
#'
#' @inheritParams load_pis
#'
#' @return Data frame, or [sf] object if `return_sf = TRUE`, containing stixel
#'   summary data. Data are organized with one stixel per row and each stixel
#'   identified by a unique `stixel_id`, the centroid of each stixel in space
#'   and time is specified by `lat`, `lon`, and `date`.
#'
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
#' dplyr::glimpse(stixels)
#' }
load_stixels <- function(path, ext, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  db_file <- file.path(path, "stixel_summary.db")
  if(!file.exists(db_file)) {
    stop("The file 'stixel_summary.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  if (missing(ext)) {
    sql <- "SELECT * FROM stixel_summary AS s;"
  } else {
    sql <- stringr::str_glue("SELECT * FROM stixel_summary AS s ",
                             "{sql_extent_subset(ext)};")
  }

  # extract from database
  stx <- DBI::dbGetQuery(db, sql)
  stx <- dplyr::tibble(stx)
  DBI::dbDisconnect(db)

  # clean names
  names(stx) <- stringr::str_to_lower(names(stx))

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(stx[, c("latitude", "longitude",
                                                "day_of_year")])
  if (any(!has_centroid)) {
    warning("Removing ", sum(!has_centroid),
            " stixels with missing centroids data.")
    stx <- stx[has_centroid, ]
  }

  # subset to polygon if one was provided
  if (!missing(ext) && ext$type == "polygon") {
    stx <- ebirdst_subset(stx, ext = ext)
  }

  if (isTRUE(return_sf)) {
    stx <- sf::st_as_sf(stx, coords = c("longitude", "latitude"), crs = 4326)
  }
  return(stx)
}


#' Load eBird Status and Trends test data predictions
#'
#' During eBird Status and Trends modeling, predictions are made for checklists
#' in a test dataset that is not included in the model fitting process. This
#' function loads these predictions in addition to the actual observed count on
#' the associated checklist. These data are used by [ebirdst_ppms()] for
#' calculating predictive performance metrics.
#'
#' @inheritParams load_pis
#'
#' @return Data frame, or [sf] object if `return_sf = TRUE`, containing
#'   observed counts and model predictions for the test data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # test data
#' test_predictions <- load_predictions(path)
#' dplyr::glimpse(test_predictions)
#' }
load_predictions <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  db_file <- file.path(path, "predictions.db")
  if(!file.exists(db_file)) {
    stop("The file 'predictions.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  preds <- DBI::dbGetQuery(db, "SELECT * FROM predictions;")
  preds <- dplyr::tibble(preds)
  DBI::dbDisconnect(db)

  if (isTRUE(return_sf)) {
    preds <- sf::st_as_sf(preds, coords = c("longitude", "latitude"),
                          crs = 4326)
  }
  return(preds)
}


#' Load eBird Status and Trends configuration file
#'
#' Load the configuration file for an eBird Status and Trends runs. This
#' configuration file is mostly for internal use and contains a variety of
#' parameters used in the modeling process.
#'
#' @inheritParams load_raster
#'
#' @return A list with the run configuration parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load configuration file
#' cfg <- load_config(path)
#' }
load_config <- function(path) {
  stopifnot(dir.exists(path))
  cfg_file <- file.path(path, "config.json")
  if(!file.exists(cfg_file)) {
    stop("The file 'config.json' does not exist in: ", path)
  }
  # load configuration file
  p <- jsonlite::read_json(cfg_file, simplifyVector = TRUE)
  names(p) <- tolower(names(p))
  return(p)
}


#' Load full annual cycle map parameters
#'
#' Get the map parameters used on the eBird Status and Trends website to
#' optimally display the full annual cycle data. This includes bins for the
#' abundance data, a projection, and an extent to map. The extent is the spatial
#' extent of non-zero data across the full annual cycle and the projection is
#' optimized for this extent.
#'
#' @inheritParams load_raster
#'
#' @return A list containing elements:
#' - `custom_projection`: a custom projection optimized for the given species'
#'    full annual cycle
#' - `fa_extent`: a [SpatExtent][terra::ext()] object storing the spatial extent of non-zero
#' data for the given species in the custom projection
#' - `res`: a numeric vector with 2 elements giving the target resolution of
#'    raster in the custom projection
#' - `fa_extent_sinu`: the extent in sinusoidal projection
#' - `weekly_bins`/`weekly_labels`: weekly abundance bins and labels for the
#' full annual cycle
#' - `seasonal_bins`/`seasonal_labels: seasonal abundance bins and labels for
#' the full annual cycle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # get map parameters
#' load_fac_map_parameters(path)
#' }
load_fac_map_parameters <- function(path) {
  stopifnot(dir.exists(path))

  # load configuration file
  p <- load_config(path)
  ext_order <- unlist(p$bbox_sinu)[c("xmin", "xmax", "ymin", "ymax")]

  seasonal_bins <- list(custom_projection = p$projection$crs,
                        fa_extent = terra::ext(p$projection$extent),
                        res = p$projection$res,
                        fa_extent_sinu = terra::ext(ext_order),
                        weekly_bins = p$bins$hr$breaks,
                        weekly_labels = p$bins$hr$labels,
                        seasonal_bins = p$bins_seasonal$hr$breaks,
                        seasonal_labels = p$bins_seasonal$hr$labels)
}


# internal ----

sql_extent_subset <- function(ext) {
  stopifnot(inherits(ext, "ebirdst_extent"))
  if (!sf::st_is_longlat(ext$extent)) {
    stop("Load functions can only subset to spatial extents defined in ",
         "unprojected, latitude-longitude coordinates. Considering loading ",
         "then subsetting with ebirdst_subset().")
  }

  # spatial filtering
  b <- sf::st_bbox(ext$extent)
  sql <- stringr::str_glue("WHERE s.longitude > {b['xmin']} ",
                           "AND s.longitude <= {b['xmax']} ",
                           "AND s.latitude > {b['ymin']} ",
                           "AND s.latitude <= {b['ymax']}")


  # temporal filtering
  if (!identical(ext$t, c(0, 1))) {
    t <- ext$t * 366
    if (t[1] <= t[2]) {
      t_sql <- stringr::str_glue("AND s.day_of_year > {t[1]} ",
                                 "AND s.day_of_year <= {t[2]}")
      sql <- paste(sql, t_sql)
    } else {
      t_sql <- stringr::str_glue("AND (s.day_of_year > {t[1]} ",
                                 "OR s.day_of_year <= {t[2]})")
      sql <- paste(sql, t_sql)
    }
  }
  return(sql)
}


transform_predictor_names <- function(x) {
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "\\.", "_")
  x[x == "i_stationary"] <- "is_stationary"
  x[x == "lon"] <- "longitude"
  x[x == "lat"] <- "latitude"
  return(x)
}
