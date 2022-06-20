#' Load eBird Status and Trends raster data cubes
#'
#' Each of the eBird Status and Trends products is packaged as a GeoTIFF file
#' (referred to as a "cube) with 52 bands, one for each week of the year. This
#' function loads the cube for a given product and species as a `RasterStack`
#' object.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#' @param product character; Status and Trends product to load, see Details for
#'   available products. It is also possible to return a template raster with no
#'   data.
#' @param resolution character; the resolution of the raster data to load. The
#'   default is to load the native ~3 km resolution (`"hr"`); however, for some
#'   applications 9 km (`"mr"`) or 27 km (`"lr"`) data may be suitable.
#'
#' @details The available Status and Trends data cubes are as follows:
#'
#' - `occurrence`: the expected probability of occurrence of the species,
#' ranging from 0 to 1, on an eBird Traveling Count by a skilled eBirder
#' starting at the optimal time of day with the optimal search duration and
#' distance that maximizes detection of that species in a region.
#' - `count`: the expected count of a species, conditional on its occurrence at
#' the given location, on an eBird Traveling Count by a skilled eBirder
#' starting at the optimal time of day with the optimal search duration and
#' distance that maximizes detection of that species in a region.
#' - `abundance`: the expected relative abundance, computed as the product of
#' the probability of occurrence and the count conditional on occurrence, of the
#' species on an eBird Traveling Count by a skilled eBirder starting at the
#' optimal time of day with the optimal search duration and distance that
#' maximizes detection of that species in a region.
#' - `abundance_lower`: the lower 10th quantile of the expected relative
#' abundance of the species on an eBird Traveling Count by a skilled eBirder
#' starting at the optimal time of day with the optimal search duration and
#' distance that maximizes detection of that species in a region.
#' - `abundance_upper`: the upper 90th quantile of the expected relative
#' abundance of the species on an eBird Traveling Count by a skilled eBirder
#' starting at the optimal time of day with the optimal search duration and
#' distance that maximizes detection of that species in a region.
#'
#' In addition to these cubes with 52 layers (one for each week), it is possible
#' to load:
#' - `abundance_seasonal`: the expected relative abundance averaged across the
#' weeks within each season. The date boundaries used for the seasonal
#' definitions appear in `ebirdst_runs` and if a season failed review no
#' associated layer will be included.
#' - `template`: a template raster covering the whole Earth and without any
#' data.
#'
#' @return A `RasterStack` with 52 layers for the given product, labeled by
#'   week. Seasonal abundance will have up to four layers labeled according to
#'   the seasons. The template raster will be returned as a `RasterLayer`.
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
#' # load data
#' load_raster(path, "abundance")
#' }
load_raster <- function(path,
                        product = c("abundance",
                                    "abundance_seasonal",
                                    "count",
                                    "occurrence",
                                    "abundance_lower",
                                    "abundance_upper",
                                    "template"),
                        resolution = c("hr", "mr", "lr")) {

  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  product <- match.arg(product)
  resolution <- match.arg(resolution)

  if (product %in% c("abundance", "count", "occurrence")) {
    product <- paste0(product, "_median")
  }

  # load raster
  if (product == "template") {
    if (resolution != "hr") {
      stop("For the raster template, resolution must be 'hr'")
    }
    # template raster
    tif_path <- file.path(path, "srd_raster_template.tif")
    if (length(tif_path) != 1 || !file.exists(tif_path)) {
      stop("Error locating the raster template")
    }
    return(suppressWarnings(raster::raster(tif_path)))
  } else if (product == "abundance_seasonal") {
    # seasonal abundance
    pat <- stringr::str_glue("abundance.*{resolution}.*\\.tif$")
    tif_path <- list.files(file.path(path, "seasonal"),
                           pattern = pat,
                           full.names = TRUE)
    # drop core seasons
    is_core <- stringr::str_detect(tif_path, "_core_")
    tif_path <- tif_path[!is_core]
    if (any(!file.exists(tif_path))) {
      stop("Error locating seasonal abundance GeoTIFFs")
    } else if (length(tif_path) == 0) {
      stop("No seasonal abundance GeoTIFFs found")
    }
    season_order <- c("breeding", "postbreeding_migration",
                      "nonbreeding", "prebreeding_migration",
                      "resident")
    seasons <- stringr::str_extract(tif_path,
                                    "(?<=_abundance_)[-a-z]+")
    seasons <- stringr::str_replace_all(seasons, "-", "_")
    r <- suppressWarnings(raster::stack(tif_path))
    names(r) <- seasons
    return(r[[intersect(season_order, seasons)]])
  } else {
    # weekly stack
    pat <- stringr::str_glue("{product}_{resolution}.*\\.tif$")
    tif_path <- list.files(file.path(path, "cubes"),
                           pattern = pat,
                           full.names = TRUE)
    if (length(tif_path) != 1 || !file.exists(tif_path)) {
      stop(paste("Error locating GeoTIFF file for:", product))
    }
    r <- suppressWarnings(raster::stack(tif_path))

    # label
    if (raster::nlayers(r) == 52) {
      r <- label_raster_stack(r)
    } else {
      l <- load_config(path)
      weeks <- paste0(l$SRD_PRED_YEAR, "-", l$DATE_NAMES)
      weeks <- as.Date(weeks, "%Y-%m-%d")
      r <- label_raster_stack(r, weeks = weeks)
    }
    return(r)
  }
}


#' Load eBird Status and Trends predictor importance data
#'
#' Loads the predictor importance (PI) data from the stixel_summary.db sqlite
#' database. PI estimates are provided for each stixel over which a model was
#' run and are identified by a unique stixel ID in addition to the coordinates
#' of the stixel centroid. PI estimates are for the occurrence model only.
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
#'   estimates for each stixel for both the occurrence and relative abundance
#'   models. The data are provided in a 'wide' format, with each row
#'   corresponding to the PI estimates for a give stixel for the occurrence
#'   count model, and the relative importance of each predictor in columns.
#'   Stixels are identified by a unique `stixel_id`, the centroid of the stixel
#'   in space and time is specified by the `lat`, `lon`, and `date` column,
#'   which expresses the day of year as a value from 0-1.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load predictor importance
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
#'   estimates for each stixel for either the occurrence and relative model. The
#'   data frame will have the following columns:
#'   - `stixel_id`: unique stixel identifier
#'   - `lat` and `lon`: stixel centroid
#'   - `date`: day of year, expressed as a value from 0-1, of the stixel center
#'   - `predictor`: name of the predictor that the PD data correspond to, for a
#'   full list of predictors consult the [ebirdst_predictors] data frame
#'   - `predictor_value`: value of the predictor variable at which PD is
#'   evaluated
#'   - `response`: predicted response, occurrence or relative abundance, at the
#'   given value of the predictor averaged across all the values of the other
#'   predictors
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load partial dependence data
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
#' \donttest{
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
#' the associated checklist. These data are used by [ebirdst_ppms()] to get for
#' calculating predictive performance metrics.
#'
#' @inheritParams load_pis
#'
#' @return Data frame, or [sf] object if `return_sf = TRUE`, containing
#'   observed counts and model predictiosn for the test data.
#'
#' @export
#'
#' @examples
#' \donttest{
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
#' \donttest{
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
    stop("The file 'config.rds' does not exist in: ", path)
  }
  # load configuration file
  jsonlite::read_json(cfg_file, simplifyVector = TRUE)
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
#' - `fa_extent`: an `Extent` object storing the spatial extent of non-zero
#'    data for the given species in the custom projection
#' - `res`: a numeric vector with 2 elements giving the target resolution of
#'    raster in the custom projection.
#' - `fa_extent_sinu`: the extent in sinusoidal projection
#' - `abundance_bins`: abundance bins for the full annual cycle
#'
#' @export
#'
#' @examples
#' \donttest{
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

  list(custom_projection = p$projection$crs,
       fa_extent = raster::extent(p$projection$extent),
       res = p$projection$res,
       fa_extent_sinu = raster::extent(unlist(p$bbox_sinu)),
       abundance_bins = p$bins$hr$breaks)
}

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
