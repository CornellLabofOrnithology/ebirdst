#' Download eBird Status and Trends Data Products
#'
#' Download an eBird Status and Trends data package for a single species, or for
#' an example species. Accessing Status and Trends data requires an access key,
#' consult [set_ebirdst_access_key()] for instructions on how to obtain and
#' store this key. The example data consist of the results for Yellow-bellied
#' Sapsucker subset to Michigan and are much smaller than the full dataset,
#' making these data quicker to download and process. Only the low resolution
#' data are available for the example data In addition, the example data are
#' accessible without an access key.
#'
#' @param species character; a single species given as a scientific name, common
#'   name or six-letter species code (e.g. woothr). The full list of valid
#'   species is can be viewed in the [ebirdst_runs] data frame included in this
#'   package. To download the example dataset, use "example_data".
#' @param path character; directory to download the data to. All downloaded
#'   files will be placed in a sub-directory of this directory named for the
#'   data version year, e.g. "2020" for the 2020 Status Data Products. Each species'
#'   data package will then appear in a directory named with the eBird species
#'   code. Defaults to a persistent data directory, which can be found by
#'   calling `ebirdst_data_dir()`.
#' @param tifs_only logical; whether to only download the GeoTIFFs for
#'   abundance and occurrence (the default), or download the entire data
#'   package, including data for predictor importance, partial dependence, and
#'   predictive performance metrics.
#' @param force logical; if the data have already been downloaded, should a
#'   fresh copy be downloaded anyway.
#' @param show_progress logical; whether to print download progress information.
#' @param pattern character; regular expression pattern to supply to
#'   [str_detect()][stringr::str_detect()] to filter files to download. Note
#'   that some files are mandatory and will always be downloaded.
#' @param dry_run logical; whether to do a dry run, just listing files that will
#'   be downloaded. This can be useful when testing the use of `pattern` to
#'   filter the files to download.
#'
#' @return Path to the folder containing the downloaded data package for the
#'   given species. If `dry_run = TRUE` a list of files to download will be
#'   returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download the example data
#' ebirdst_download("example_data")
#'
#' # download the data package for wood thrush, geotiffs only
#' ebirdst_download("woothr")
#' # download the data package for wood thrush, all data
#' ebirdst_download("woothr", tifs_only = FALSE)
#'
#' # use pattern to only download low resolution data
#' # dry_run can be used to see what files will be downloaded
#' ebirdst_download("lobcur", pattern = "_lr_", dry_run = TRUE)
#' # use pattern to only download the high res weekly abundance data
#' ebirdst_download("lobcur", pattern = "abundance_median_hr", dry_run = TRUE)
#' }
ebirdst_download <- function(species,
                             path = ebirdst_data_dir(),
                             tifs_only = TRUE,
                             force = FALSE,
                             show_progress = TRUE,
                             pattern = NULL,
                             dry_run = FALSE) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.logical(tifs_only), length(tifs_only) == 1)
  stopifnot(is.logical(force), length(force) == 1)
  stopifnot(is.logical(show_progress), length(show_progress) == 1)
  stopifnot(is.logical(dry_run), length(dry_run) == 1)

  species <- tolower(species)
  is_example <- (species == "example_data")

  # example data or a real run
  if (is_example) {
    api_url <- paste0("https://raw.githubusercontent.com/",
                     "ebird/ebirdst_example-data/main/",
                     "example-data/")
    # file list
    fl <- system.file("extdata", "example-data_file-list.txt",
                      package = "ebirdst")
    files <- readLines(fl)
    year_species <- stringr::str_split(files[1], "/")[[1]][1:2]
    version_year <- year_species[1]
    species <- year_species[2]
  } else {
    # version of the data products that this package version corresponds to
    version_year <- ebirdst_version()[["version_year"]]
    species <- get_species(species)
    if (is.na(species)) {
      stop("The requested species was not modeled by Status and Trends. ",
           "Consult ebirdst_runs for a complete list of available species.")
    }

    # api url and key
    key <- get_ebirdst_access_key()
    api_url <- "https://st-download.ebird.org/v1"

    # get file list for this species
    list_obj_url <- stringr::str_glue("{api_url}/list-obj/{version_year}/",
                                      "{species}?key={key}")
    files <- tryCatch(suppressWarnings({
      jsonlite::read_json(list_obj_url, simplifyVector = TRUE)
    }), error = function(e) NULL)
    if (is.null(files)) {
      stop("Cannot access Status and Trends data URL. Ensure that you have ",
           "a working internet connection and a valid API key for the Status ",
           "and Trends data. Note that the API keys expire after 1 month, so ",
           "may need to update your key. Visit https://ebird.org/st/request")
    }

    # remove web_download folder
    web_down <- stringr::str_detect(dirname(files), pattern = "web_download")
    files <- files[!web_down]
  }

  if (length(files) == 0) {
    stop("No data found for species ", species)
  }

  # path to data package
  run_path <- file.path(path, version_year, species)

  # only download databases when explicitly requested
  if (isTRUE(tifs_only)) {
    files <- files[!stringr::str_detect(files, "\\.db$")]
  }

  # apply pattern
  if (!is.null(pattern)) {
    stopifnot(is.character(pattern), length(pattern) == 1, !is.na(pattern))
    pat_match <- stringr::str_detect(basename(files), pattern = pattern)
    if (all(!pat_match)) {
      stop("No files matched pattern")
    }

    # always download config file
    is_config <- stringr::str_detect(basename(files), pattern = "config.json$")
    files <- files[pat_match | is_config]
  }

  # print files to download for dry run
  if (dry_run) {
    message("Downloading data package for ", species, " to:\n  ", path)
    message(paste(c("File list:", files), collapse = "\n  "))
    return(invisible(files))
  }

  # prepare download paths
  files <- data.frame(file = files)
  if (is_example) {
    files$src_path <- paste0(api_url, files$file)
  } else {
    files$src_path <- stringr::str_glue("{api_url}/fetch?objKey={files$file}",
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
      return(invisible(normalizePath(run_path)))
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
  for (i in seq_len(nrow(files))) {
    dl_response <- utils::download.file(files$src_path[i],
                                        files$dest_path[i],
                                        quiet = !show_progress,
                                        mode = "wb")
    if (dl_response != 0) {
      stop("Error downloading file: ", files$file[i])
    }
  }
  options(timeout = old_timeout)
  return(invisible(normalizePath(run_path)))
}


#' Get the path to the data package for a given species
#'
#' This helper function can be used to get the path to a data package for a
#' given species to be used by the various data loading functions.
#'
#' @inheritParams ebirdst_download
#'
#' @return The path to the data package directory.
#' @export
#'
#' @examples
#' \dontrun{
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
#' path <- get_species_path("Yellow-bellied Sapsucker")
#' path <- get_species_path("Sphyrapicus varius")
#' path <- get_species_path("yebsap")
#' }
get_species_path <- function(species, path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))

  # append version to path
  path <- file.path(path, ebirdst_version()["version_year"])

  if (species == "example_data") {
    run <- "yebsap-example"
  } else {
    species <- get_species(species)
    row_id <- which(ebirdst::ebirdst_runs$species_code == species)
    if (length(row_id) != 1) {
      stop(sprintf("species = %s does not uniquely identify a species.",
                   species))
    }
    run <- ebirdst::ebirdst_runs$species_code[row_id]
  }
  species_path <- path.expand(file.path(path, run))
  if (!dir.exists(species_path)) {
    stop(paste("No data package found for species:", species))
  }
  return(species_path)
}


#' Path to eBird Status and Trends data download directory
#'
#' Identify and return the path to the default download directory for eBird
#' Status and Trends data products. This directory can be defined by setting the
#' environment variable `EBIRDST_DATA_DIR`, otherwise the directory returned by
#' `tools::R_user_dir("ebirdst", which = "data")` will be used.
#'
#' @return The path to the data download directory.
#' @export
#'
#' @examples
#' ebirdst_data_dir()
ebirdst_data_dir <- function() {
  env_var <- Sys.getenv("EBIRDST_DATA_DIR")
  if (is.null(env_var) || env_var == "") {
    return(tools::R_user_dir("ebirdst", which = "data"))
  } else {
    return(env_var)
  }
}


#' eBird Status and Trends Data Products version
#'
#' Identify the version of the eBird Status and Trends Data Products that this
#' version of the R package works with. Versions are defined by the year that
#' all model estimates are made for. In addition, the release data and end date
#' for access of this version of the data are provided. Note that after the
#' given access end data you will no longer be able to download this version of
#' the data and will be required to update the R package and transition to using
#' a newer data version.
#'
#' @return A list with three components: `version_year` is the year the model
#'   estimates are made for in this version of the data, `release_year` is the
#'   year this version of the data were released, and `access_end_date` is the
#'   last date that users will be able to download this version of the data.
#' @export
#'
#' @examples
#' ebirdst_version()
ebirdst_version <- function() {
  list(version_year = 2021,
       release_year = 2022,
       access_end_date = as.Date("2023-11-30"))
}
