#' Download eBird Status and Trends Data
#'
#' Download an eBird Status and Trends data package for a single species, or for
#' an example species, to a specified path. The example data consist of the
#' results for Yellow-bellied Sapsucker subset to Michigan and are much smaller
#' than the full dataset, making these data quicker to download and process.
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
#'
#' @return Path to the run-specific root of the downloaded files.
#'
#' @export
#'
#' @examples
#' # download the example data
#' ebirdst_download("example_data")
#'
#' \dontrun{
#' # download the full data package for wood thrush
#' ebirdst_download("woothr")
#' }
ebirdst_download <- function(species,
                             path = rappdirs::user_data_dir("ebirdst"),
                             tifs_only = TRUE,
                             force = FALSE) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.logical(tifs_only), length(tifs_only) == 1)
  stopifnot(is.logical(force), length(force) == 1)
  species <- tolower(species)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # example data or a real run
  bucket_url <- "https://ebirdst-data.s3.amazonaws.com/"
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
  bucket_url_sp <- paste0(bucket_url, "?prefix=", run)

  # get bucket contents
  s3_contents <- xml2::xml_ns_strip(xml2::read_xml(bucket_url_sp))
  s3_contents <- xml2::xml_find_all(s3_contents, ".//Contents")
  if (length(s3_contents) == 0) {
    stop(sprintf("Files not found on AWS S3 for species = %s", species))
  }

  # store filename and size
  s3_files <- data.frame(
    file = xml2::xml_text(xml2::xml_find_all(s3_contents, ".//Key")),
    size = xml2::xml_text(xml2::xml_find_all(s3_contents, ".//Size")),
    stringsAsFactors = FALSE)
  s3_files$size <- as.numeric(s3_files$size)

  # filter to desired run/species
  s3_files <- s3_files[as.numeric(s3_files$size) > 0 &
                         grepl(run, s3_files$file), ]
  if (nrow(s3_files) == 0) {
    stop(sprintf("Files not found on AWS S3 for species = %s", species))
  }

  # only dl rasters unless requested otherwise
  if (isTRUE(tifs_only)) {
    s3_files <- s3_files[grepl("tif$", s3_files$file) |
                           grepl("rds$", s3_files$file) |
                           grepl("band_dates", s3_files$file), ]
  }

  dl_filter <- Sys.getenv("EBIRDST_DL_FILTER")
  if (dl_filter != "") {
    s3_files <- s3_files[grepl(dl_filter, s3_files$file), ]
  }

  # human readable download size if we want to add a message
  #size_human <- utils:::format.object_size(sum(s3_files$size), "auto")

  # prepare download paths
  s3_files$s3_path <- paste0(bucket_url, s3_files$file)
  s3_files$local_path <- file.path(path, s3_files$file)
  s3_files$exists <- file.exists(s3_files$local_path)
  # create necessary directories
  dirs <- unique(dirname(s3_files$local_path))
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # check if already exists
  if (all(s3_files$exists)) {
    if (!force) {
      message("Data already exists, use force = TRUE to re-download.")
      return(invisible(normalizePath(file.path(path, run))))
    }
  } else if (any(s3_files$exists)) {
    if (!force) {
      message(paste("Some files already exist, only downloading new files.",
                    "\nUse force = TRUE to re-download all files."))
      s3_files <- s3_files[!s3_files$exists, ]
    }
  }

  # download
  for(i in seq_len(nrow(s3_files))) {
    dl_response <- utils::download.file(s3_files$s3_path[i],
                                        s3_files$local_path[i],
                                        quiet = TRUE,
                                        mode = "wb")
    if (dl_response != 0) {
      stop("Error downloading files from AWS S3")
    }
  }

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
#' sp_path <- get_species_path("example_data")
#'
#' # use it to load data
#' abd <- load_raster("abundance", sp_path)
#'
#' # get the path to the full data package for yellow-bellied sapsucker
#' # common name, scientific name, or species code can be used
#' sp_path <- get_species_path("Yellow-bellied Sapsucker")
#' sp_path <- get_species_path("Sphyrapicus varius")
#' sp_path <- get_species_path("yebsap")
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
  species_path <- file.path(path, run)
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
#' sp_path <- ebirdst_download("example_data")
#' # or get the path if you already have the data downloaded
#' sp_path <- get_species_path("example_data")
#'
#' # load data
#' load_raster(sp_path, "abundance")
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
    tif_path <- list.files(file.path(path, "data"),
                           pattern = "srd_raster_template\\.tif$",
                           full.names = TRUE)
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
#' centroid.
#'
#' @inheritParams load_raster
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data frame containing PI estimates for each stixel, as well as stixel
#' summary information.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' sp_path <- get_species_path("example_data")
#'
#' # load predictor importance
#' pis <- load_pis(sp_path)
#'
#' # plot the top 15 predictor importances
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pis(pis, ext = e, n_top_pred = 15, by_cover_class = TRUE)
#' }
load_pis <- function(path, type = c("occurrence", "abundance"),
                     return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  db_file <- file.path(path, "pi-pd.db")
  if(!file.exists(db_file)) {
    stop("The file 'pi-pd.db' does not exist in: ", path)
  }

  # load pi file
  pi_df <- data.table::fread(pi_file,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  pi_df <- pi_df[, c("stixel", "data_type") := NULL]

  # summary file
  summary_df <- load_summary(path)
  summary_df <- summary_df[, 1:12]

  # merge pis with summary
  pi_summary <- dplyr::inner_join(summary_df, pi_df, by = "stixel_id")
  pi_summary <- as.data.frame(pi_summary)
  names(pi_summary) <- tolower(names(pi_summary))

  if (isTRUE(return_sf)) {
    pi_summary <- sf::st_as_sf(pi_summary, coords = c("lon", "lat"), crs = 4326)
  }
  return(pi_summary)
}


#' Load partial dependency data for the Status and Trends predictors
#'
#' Partial dependence (PD) plots depict the relationship between the modeled
#' response (occurrence or relative abundance) and the predictor variables used
#' in the model. Status and Trends provides the data to generate these plots for
#' every stixel.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing predictor importance values for each stixel, as
#'   well as stixel summary information.
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # load predictor dependencies
#' pds <- load_pds(sp_path)
#'
#' # plot the top 15 predictor importances
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pds(pds, ext = e, n_top_pred = 15, by_cover_class = TRUE)
load_pds <- function(path, model = c("occurrence", "abundance"),
                     return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)
  stopifnot(model %in% c("occurrence", "abundance"))

  stixel_path <- file.path(path, "results", "stixels")
  pd_file <- file.path(stixel_path, "pd.txt")

  if(!file.exists(pd_file)) {
    stop(paste("The file pi.txt does not exist at:", stixel_path))
  }

  # load pi file
  pd_df <- data.table::fread(pd_file,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)

  # pick occ or abd pis
  model <- match.arg(model)
  if (model == "occurrence") {
    pd_df <- pd_df[pd_df$model == "occ", ]
  } else {
    pd_df <- pd_df[pd_df$model == "abd", ]
  }

  # summary file
  summary_df <- load_summary(path)
  summary_df <- summary_df[, 1:5]

  # merge pis with summary
  pd_summary <- dplyr::inner_join(summary_df, pd_df, by = "stixel_id")
  pd_summary <- as.data.frame(pd_summary)
  names(pd_summary) <- tolower(names(pd_summary))

  if (isTRUE(return_sf)) {
    pd_summary <- sf::st_as_sf(pd_summary, coords = c("lon", "lat"), crs = 4326)
  }
  return(pd_summary)
}


#' Test data loader
#'
#' Internal function used by [compute_ppms()] to get the full test data set for
#' calculating predictive performance metrics. This file contains the observed
#' counts and all predictor variables used for modeling for each checklist in
#' the test dataset.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing test data, including checklist locations,
#'   observed counts, and predictor variables.
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # test data
#' test_data <- ebirdst:::load_test_data(sp_path)
#' dplyr::glimpse(test_data)
#' }
load_test_data <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  data_path <- file.path(path, "data")
  td_file <- list.files(file.path(path, "data"),
                        pattern = "test-data\\.csv$",
                        full.names = TRUE)
  if (length(td_file) != 1 || !file.exists(td_file)) {
    stop(paste("Error locating raw test data file at:", data_path))
  }

  # load raw test data from sqlite db
  td_df <- data.table::fread(td_file,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  td_df[["DATA_TYPE"]] <- NULL

  # fix names
  names(td_df) <- tolower(names(td_df))
  nm_idx <- match(c("longitude", "latitude"), names(td_df))
  names(td_df)[nm_idx] <- c("lon", "lat")

  td_df <- as.data.frame(td_df)
  if (isTRUE(return_sf)) {
    td_df <- sf::st_as_sf(td_df, coords = c("lon", "lat"), crs = 4326)
  }
  return(td_df)
}


#' Test data predictions loader
#'
#' Loads the model predictions for each checklist in the test dataset. Median,
#' and upper and lower confidence intervals are provided for predicted
#' occurrence, count, and relative abundance.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing median, and upper and lower confidence
#'   intervals are provided for predicted occurrence, count, and relativce
#'   abundance.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # test data
#' test_preds <- load_test_preds(sp_path)
#' dplyr::glimpse(test_preds)
#' }
load_test_preds <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  pred_path <- file.path(path, "results", "preds")
  td_file <- file.path(pred_path, "test_pred_ave.txt")

  if(!file.exists(td_file)) {
    stop(paste("The file test_pred_ave.txt does not exist at:", pred_path))
  }

  # load test data
  td_df <- data.table::fread(td_file,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  td_df <- td_df[, "data_type" := NULL]
  td_df <- as.data.frame(td_df, stringsAsFactors = FALSE)

  if (isTRUE(return_sf)) {
    td_df <- sf::st_as_sf(td_df, coords = c("lon", "lat"), crs = 4326)
  }

  return(td_df)
}


#' Load configuration file
#'
#' Internal function used by load_summary(), load_pis(), and load_pds() to get
#' configuration variables from STEM species run information
#' (from *_config.RData).
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#'
#' @return environment object containing all run parameters.
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' ebirdst:::load_config(sp_path)
#' }
load_config <- function(path) {
  config_file <- list.files(paste(path, "/data", sep = ""),
                            pattern = "*_config.rds",
                            full.names = TRUE)

  if(length(config_file) != 1 || !file.exists(config_file)) {
    stop(paste("*_config.rds file does not exist in the /data directory. ",
               "Check your paths so that they look like this: ",
               "~/directory/<six_letter_code-ERD2016-PROD-date-uuid>/. ",
               "Make sure you do not change the file structure of the results.",
               sep = ""))
  }

  return(readRDS(config_file))
}


#' Load full annual cycle map parameters
#'
#' Get the map parameters used on the eBird Status and Trends website to
#' optimally display the full annual cycle data. This includes bins for the
#' abundance data, a projection, and an extent to map. The extent is the spatial
#' extent of non-zero data across the full annual cycle and the projection is
#' optimized for this extent.
#'
#' @param path character; directory that the Status and Trends data for a given
#'   species was downloaded to. This path is returned by `ebirdst_download()`
#'   or `get_species_path()`.
#'
#' @return A list containing elements:
#' - `custom_projection`: a custom projection optimized for the given species'
#'    full annual cycle
#' - `fa_extent`: an `Extent` object storing the spatial extent of non-zero
#'    data for the given species in the custom projection
#' - `fa_extent_sinu`: the extent in sinusoidal projection
#' - `abundance_bins`: abundance bins for the full annual cycle
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # get map parameters
#' load_fac_map_parameters(sp_path)
#' }
load_fac_map_parameters <- function(path) {
  l <- load_config(path)

  list(custom_projection = l$CUS_PRJ,
       fa_extent = l$FA_EXTENT,
       fa_extent_sinu = l$SINU_EXTENT,
       abundance_bins = l$ABUND_BINS)
}


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
