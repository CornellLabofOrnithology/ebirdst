#' Download eBird Status and Trends Data
#'
#' Download an eBird Status and Trends data package for a single species, or for
#' an example species, to a specified path. Accessing Status and Trends data
#' requires an access key, consult [set_ebirdst_access_key()] for instructions
#' on how to obtain and store this key. The example data consist of the results
#' for Yellow-bellied Sapsucker subset to Michigan and are much smaller than the
#' full dataset, making these data quicker to download and process. In addition,
#' the example data are accessible without an access key.
#'
#' @param species character; a single species given as a scientific name, common
#'   name or six-letter species code (e.g. woothr). The full list of valid
#'   species is can be viewed in the [ebirdst_runs] data frame included in this
#'   package. To download the example dataset, use "example_data".
#' @param path character; directory to download the data to. All downloaded
#'   files will be placed in a sub-directory of this directory named according
#'   to the unique run ID associated with this species. Defaults to a persistent
#'   data directory, which can be found by calling
#'   rappdirs::user_data_dir("ebirdst")).
#' @param tifs_only logical; whether to only download the GeoTIFFs for
#'   abundance and occurrence (the default), or download the entire data
#'   package, including data for predictor importance, partial dependence, and
#'   predictive performance metrics.
#' @param force logical; if the data have already been downloaded, should a
#'   fresh copy be downloaded anyway.
#' @param show_progress logical; whether to print download progress information.
#'
#' @return Path to the folder containing the downloaded data package for the
#'   given species.
#'
#' @export
#'
#' @examples
#' # download the example data
#' ebirdst_download("example_data")
#'
#' \dontrun{
#' # download the data package for wood thrush, geotiffs only
#' ebirdst_download("woothr")
#' # download the data package for wood thrush, all data
#' ebirdst_download("woothr", tifs_only = FALSE)
#' }
ebirdst_download <- function(species,
                             path = rappdirs::user_data_dir("ebirdst"),
                             tifs_only = TRUE,
                             force = FALSE,
                             show_progress = TRUE) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.logical(tifs_only), length(tifs_only) == 1)
  stopifnot(is.logical(force), length(force) == 1)
  stopifnot(is.logical(show_progress), length(show_progress) == 1)
  species <- tolower(species)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # example data or a real run
  if (species == "example_data") {
    api_url <- "https://ebirdst-data.s3.amazonaws.com/"
    run <- "yebsap-ERD2019-STATUS-20200930-8d36d265-example"
    bucket_url_sp <- paste0(api_url, "?prefix=", run)

    # get bucket contents
    s3_contents <- xml2::xml_ns_strip(xml2::read_xml(bucket_url_sp))
    s3_contents <- xml2::xml_find_all(s3_contents, ".//Contents")
    if (length(s3_contents) == 0) {
      stop(sprintf("Files not found on AWS S3 for species = %s", species))
    }

    # store filename and size
    files <- data.frame(
      file = xml2::xml_text(xml2::xml_find_all(s3_contents, ".//Key")),
      size = xml2::xml_text(xml2::xml_find_all(s3_contents, ".//Size")))
    files$size <- as.numeric(files$size)

    # filter to desired run/species
    files <- files[as.numeric(files$size) > 0 & grepl(run, files$file), ,
                   drop = FALSE]
  } else {
    species <- get_species(species)
    which_run <- which(ebirdst::ebirdst_runs$species_code == species)
    run <- ebirdst::ebirdst_runs$run_name[which_run]
    if (is.na(run) || length(run) != 1) {
      stop("species does not uniquely identify a Status and Trends run.")
    }

    # api url and key
    key <- get_ebirdst_access_key()
    api_url <- "http://3.87.12.94:8080/v1/"

    # get file list for this species
    list_obj_url <- stringr::str_glue("{api_url}/list-obj/{species}?key={key}")
    files <- tryCatch(suppressWarnings({
      jsonlite::read_json(list_obj_url, simplifyVector = TRUE)
    }), error = function(e) NULL)
    if (is.null(files)) {
      stop("Cannot access Status and Trends data URL. Ensure that you have ",
           "a working internet connection and a valid API key for the Status ",
           "and Trends data.")
    }
    files <- data.frame(file = files)
  }
  if (nrow(files) == 0) {
    stop("No data found for species ", species)
  }

  # only download databases when explicitly requested
  if (isTRUE(tifs_only)) {
    files <- files[!stringr::str_detect(files$file, "\\.db$"), , drop = FALSE]
  }

  # prepare download paths
  if (species == "example_data") {
    files$src_path <- paste0(api_url, files$file)
  } else {
    files$src_path <- stringr::str_glue("{api_url}fetch?objKey={files$file}",
                                        "&key={key}")
  }
  files$dest_path <- file.path(path, files$file)
  files$exists <- file.exists(files$dest_path)
  # create necessary directories
  dirs <- unique(dirname(files$dest_path))
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # check if already exists
  if (all(files$exists)) {
    if (!isTRUE(force)) {
      message("Data already exists, use force = TRUE to re-download.")
      return(invisible(normalizePath(file.path(path, run))))
    }
  } else if (any(files$exists)) {
    if (!isTRUE(force)) {
      message(paste("Some files already exist, only downloading new files.",
                    "\nUse force = TRUE to re-download all files."))
      files <- files[!files$exists, ]
    }
  }

  # download
  old_timeout <- getOption("timeout")
  options(timeout = max(3000, old_timeout))
  for(i in seq_len(nrow(files))) {
    dl_response <- utils::download.file(files$src_path[i],
                                        files$dest_path[i],
                                        quiet = !show_progress,
                                        mode = "wb")
    if (dl_response != 0) {
      stop("Error downloading file: ", files$file[i])
    }
  }
  options(timeout = old_timeout)
  return(invisible(normalizePath(file.path(path, run))))
}


#' Get the data package path for a given species
#'
#' This helper function can be used to get the path to a data package for a
#' given species to be used by the various loading functions.
#'
#' @inheritParams ebirdst_download
#'
#' @return The path to the data package directory.
#' @export
#'
#' @examples
#' # download the example data
#' ebirdst_download("example_data")
#'
#' # get the path
#' path <- get_species_path("example_data")
#'
#' # use it to load data
#' abd <- load_raster(path, "abundance")
#'
#' # get the path to the full data package for yellow-bellied sapsucker
#' # common name, scientific name, or species code can be used
#' \dontrun{
#' path <- get_species_path("Yellow-bellied Sapsucker")
#' path <- get_species_path("Sphyrapicus varius")
#' path <- get_species_path("yebsap")
#' }
get_species_path <- function(species,
                             path = rappdirs::user_data_dir("ebirdst")) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))

  if (species == "example_data") {
    run <- "yebsap-ERD2019-STATUS-20200930-8d36d265-example"
  } else {
    species <- get_species(species)
    row_id <- which(ebirdst::ebirdst_runs$species_code == species)
    if (length(row_id) != 1) {
      stop(sprintf("species = %s does not uniquely identify a species.",
                   species))
    }
    run <- ebirdst::ebirdst_runs$run_name[row_id]
  }
  species_path <- path.expand(file.path(path, run))
  if (!dir.exists(species_path)) {
    stop(paste("No data package found for species:", species))
  }
  return(species_path)
}


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
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load data
#' load_raster(path, "abundance")
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
    tif_path <- list.files(file.path(path, "abundance_seasonal"),
                           pattern = paste0("_", resolution, "_",
                                            ".*_abundance-seasonal_.*tif$"),
                           full.names = TRUE)
    if (any(!file.exists(tif_path))) {
      stop("Error locating seasonal abundance GeoTIFFs")
    } else if (length(tif_path) == 0) {
      stop("No seasonal abundance GeoTIFFs found")
    }
    season_order <- c("breeding", "postbreeding_migration",
                      "nonbreeding", "prebreeding_migration",
                      "resident")
    seasons <- stringr::str_extract(tif_path,
                                    "(?<=abundance-seasonal_)[a-z_]+")
    r <- suppressWarnings(raster::stack(tif_path))
    names(r) <- seasons
    return(r[[intersect(season_order, seasons)]])
  } else {
    # 52 week stack
    tif_path <- list.files(file.path(path, "weekly_cubes"),
                           pattern = paste0("_", resolution, "_",
                                            ".*", product, "\\.tif$"),
                           full.names = TRUE)
    if (length(tif_path) != 1 || !file.exists(tif_path)) {
      stop(paste("Error locating GeoTIFF file for:", product))
    }
    r <- suppressWarnings(raster::stack(tif_path))
    return(label_raster_stack(r))
  }
}


#' Load eBird Status and Trends predictor importance data
#'
#' Loads the predictor importance (PI) data from the pi-pd.db sqlite database.
#' PI estimates are provided for each stixel over which a model was run and are
#' identified by a unique stixel ID in addition to the coordinates of the stixel
#' centroid. PI estimates are for the occurrence model only.
#'
#' @inheritParams load_raster
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to filter the
#'   data to. The spatial component of the extent object must be provided in
#'   unprojected, latitude-longitude coordinates.
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
load_pis <- function(path, ext, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  db_file <- file.path(path, "pi-pd.db")
  if(!file.exists(db_file)) {
    stop("The file 'pi-pd.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  if (missing(ext)) {
    sql <- paste("SELECT p.*, s.lat, s.lon, s.date",
                 "FROM occ_pis AS p INNER JOIN stixel_summary AS s",
                 "ON p.stixel_id = s.stixel_id;")
  } else {
    sql <- stringr::str_glue("SELECT p.*, s.lat, s.lon, s.date ",
                             "FROM occ_pis AS p ",
                             "INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id ",
                             "{sql_extent_subset(ext)};")
  }

  # extract from database
  pis <- DBI::dbGetQuery(db, sql)
  pis <- dplyr::tibble(pis)
  DBI::dbDisconnect(db)

  # clean up
  p <- ebirdst::ebirdst_predictors
  preds <- intersect(names(pis), p$predictor)
  pis <- pis[, c("stixel_id", "lat", "lon", "date", preds)]
  pis <- dplyr::tibble(pis)
  matches <- match(names(pis), p$predictor)
  for (i in seq_along(matches)) {
    if (!is.na(matches[i])) {
      names(pis)[i] <- p$predictor_tidy[matches[i]]
    }
  }

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(pis[, c("lat", "lon", "date")])
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
    pis <- sf::st_as_sf(pis, coords = c("lon", "lat"), crs = 4326)
  }

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
#' plot_pds(pds, "solar_noon_diff", ext = e, n_bs = 5)
#' }
load_pds <- function(path, ext, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  db_file <- file.path(path, "pi-pd.db")
  if(!file.exists(db_file)) {
    stop("The file 'pi-pd.db' does not exist in: ", path)
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(), db_file)

  # query
  if (missing(ext)) {
    sql <- paste("SELECT p.*, s.lat, s.lon, s.date",
                 "FROM occ_pds AS p INNER JOIN stixel_summary AS s",
                 "ON p.stixel_id = s.stixel_id;")
  } else {
    sql <- stringr::str_glue("SELECT p.*, s.lat, s.lon, s.date ",
                             "FROM occ_pds AS p ",
                             "INNER JOIN stixel_summary AS s ",
                             "ON p.stixel_id = s.stixel_id ",
                             "{sql_extent_subset(ext)};")
  }

  # extract from database
  pds <- DBI::dbGetQuery(db, sql)
  pds <- dplyr::tibble(pds)
  DBI::dbDisconnect(db)

  # clean up names
  preds <- ebirdst::ebirdst_predictors[, c("predictor", "predictor_tidy")]
  names(preds) <- c("covariate", "predictor")
  pds <- dplyr::inner_join(pds, preds, by = "covariate")
  pds <- pds[, c("stixel_id", "lat", "lon", "date",
                 "predictor", "predictor_value", "response")]

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(pds[, c("lat", "lon", "date")])
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
    pds <- sf::st_as_sf(pds, coords = c("lon", "lat"), crs = 4326)
  }

  return(pds)
}


#' Load summary data for eBird Status and Trends stixels
#'
#' eBird Status and Trends divides space and time into variably sized "stixels"
#' within which individual base models are fit. The process of stixelization is
#' performed many times and the prediction at any given point is the median of
#' the predictions from all the stixels that that point falls in. This function
#' loads summary statistics for each stixel, for example, the size of the
#' stixels and the number of observations within each stixel.
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

  db_file <- file.path(path, "pi-pd.db")
  if(!file.exists(db_file)) {
    stop("The file 'pi-pd.db' does not exist in: ", path)
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

  # clean up
  for (col in c("summary_hash", "stixel", "data_type", "srd_n")) {
    stx[[col]] <- NULL
  }

  # check for missing stixels centroid
  has_centroid <- stats::complete.cases(stx[, c("lat", "lon", "date")])
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
    stx <- sf::st_as_sf(stx, coords = c("lon", "lat"), crs = 4326)
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

  # strip year from date
  preds$date <- preds$date %% 1

  # clean up
  names(preds)[names(preds) == "latitude"] <- "lat"
  names(preds)[names(preds) == "longitude"] <- "lon"

  if (isTRUE(return_sf)) {
    preds <- sf::st_as_sf(preds, coords = c("lon", "lat"), crs = 4326)
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
#' # download example data
#' path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # load configuration file
#' cfg <- load_config(path)
load_config <- function(path) {
  stopifnot(dir.exists(path))
  cfg_file <- file.path(path, "config.rds")
  if(!file.exists(cfg_file)) {
    stop("The file 'config.rds' does not exist in: ", path)
  }
  # load configuration file
  readRDS(cfg_file)
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
  l <- load_config(path)

  list(custom_projection = l$FA_PROJECTION$crs,
       fa_extent = l$FA_PROJECTION$extent,
       res = l$FA_PROJECTION$res,
       fa_extent_sinu = l$SINU_EXTENT,
       abundance_bins = l$bins$hr$quantile)
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
  sql <- stringr::str_glue("WHERE ",
                           "s.lon > {b['xmin']} AND s.lon <= {b['xmax']} AND ",
                           "s.lat > {b['ymin']} AND s.lat <= {b['ymax']}")


  # temporal filtering
  t <- ext$t
  if (!identical(t, c(0, 1))) {
    if (t[1] <= t[2]) {
      t_sql <- stringr::str_glue("AND s.date > {t[1]} AND s.date <= {t[2]}")
      sql <- paste(sql, t_sql)
    } else {
      t_sql <- stringr::str_glue("AND (s.date > {t[1]} OR s.date <= {t[2]})")
      sql <- paste(sql, t_sql)
    }
  }
  return(sql)
}
