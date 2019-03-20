#' Download eBird Status and Trends Data
#'
#' Download an eBird Status and Trends data package for a single species, or for
#' an example species, to a specified path. The example data consist of the
#' results for Yellow-bellied Sapsucker subset to Michigan and are much smaller
#' than the full dataset making these data quicker to download and process.
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
#' @export
#' @examples
#' \dontrun{
#'
#' # download the example data
#' download_data("example_data")
#'
#' # download the full data package for wood thrush
#' download_data("woothr")
#'
#' }
download_data <- function(species,
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
  if (species == "example_data") {
    bucket_url <- "https://clo-is-da-example-data.s3.amazonaws.com/"
    run <- "yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83"
  } else {
    species <- get_species(species)
    row_id <- which(ebirdst::ebirdst_runs$species_code == species)
    if (length(row_id) != 1) {
      stop(sprintf("species = %s does not uniquely identify a species.",
                   species))
    }
    # presumably this url will change
    bucket_url <- "https://ebirdst-data.s3.amazonaws.com/"
    run <- ebirdst::ebirdst_runs$run_name[row_id]
  }

  if (dir.exists(file.path(path, run))) {
    if (!force) {
      message("Data already exists, use force = TRUE to re-download.")
      return(invisible(normalizePath(file.path(path, run))))
    }
  }

  # get bucket contents
  s3_contents <- xml2::xml_ns_strip(xml2::read_xml(bucket_url))
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
  s3_files[as.numeric(s3_files$size) > 0 & grepl(run, s3_files$file), ]
  if (nrow(s3_files) == 0) {
    stop(sprintf("Files not found on AWS S3 for species = %s", species))
  }

  # only dl rasters unless requested otherwise
  if (isTRUE(tifs_only)) {
    s3_files <- s3_files[grepl("results/tifs", s3_files$file) |
                           grepl("tif$", s3_files$file), ]
  }

  # human readable download size if we want to add a message
  #size_human <- utils:::format.object_size(sum(s3_files$size), "auto")

  # prepare download paths
  s3_files$s3_path <- paste0(bucket_url, s3_files$file)
  s3_files$local_path <- file.path(path, s3_files$file)
  # create necessary directories
  dirs <- unique(dirname(s3_files$local_path))
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  # download
  for(f in 1:nrow(s3_files)) {
    dl_response <- utils::download.file(s3_files[f, ]$s3_path,
                                        s3_files[f, ]$local_path,
                                        quiet = TRUE)
    if (dl_response != 0) {
      stop("Error downloading files from AWS S3")
    }
  }

  return(invisible(normalizePath(file.path(path, run))))
}


#' Load eBird Status and Trends raster data
#'
#' Each of the eBird Status and Trends products is packaged as a GeoTIFF file
#' with 52 bands, one for each week of the year. This function loads the data
#' for a given product and species as a `RasterStack` object.
#'
#' @param product character; status and trends product to load, options are
#'   relative abundance, occurrence, and upper and lower bounds on relative
#'   abundance. It is also possible to return a template raster with no data.
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#'
#' @return A `RasterStack` of data for the given product.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # load data
#' load_raster("abundance_umean", sp_path)
#'
#' }
load_raster <- function(product = c("abundance_umean",
                                    "occurrence_umean",
                                    "abundance_lower",
                                    "abundance_upper",
                                    "template"),
                        path) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  product <- match.arg(product)

  # find the file
  if (product == "template") {
    tif_path <- list.files(file.path(path, "data"),
                           pattern = "srd_raster_template\\.tif$",
                           full.names = TRUE)
  } else {
    tif_path <- list.files(file.path(path, "results", "tifs"),
                           pattern = paste0("hr_2016_", product, "\\.tif$"),
                           full.names = TRUE)
  }
  if (length(tif_path) != 1 || !file.exists(tif_path)) {
    stop(paste("Error locating GeoTIFF file for:", product))
  }

  # return
  if (product == "template") {
    return(raster::raster(tif_path))
  } else {
    return(label_raster_stack(raster::stack(tif_path)))
  }
}


#' Labels 52 week RasterStack with the dates for each band
#'
#' The `raster` package does not allow layer names to be saved with the bands of
#' a multi-band GeoTiff. Accordingly, all eBird Status and Trends products
#' raster results cover the entire 52 week temporal extent of analysis. For
#' convenience, this function labels the RasterStack once it has been loaded
#' with the dates for each band.
#'
#' @param x `RasterStack` or `RasterBrick`; original eBird Status and Trends
#'   product raster GeoTiff with 52 bands, one for each week.
#'
#' @return A `RasterStack` or `RasterBrick` with names assigned for the dates in
#'   the format of "XYYYY.MM.DD" per raster package constraints. The Raster*
#'   objects do not allow the names to start with a number, nor are they allowed
#'   to contain "-", so it is not possible to store the date in an ISO compliant
#'   format.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example abundance data
#' sp_path <- download_data("example_data")
#' abd <- load_raster("abundance_umean", sp_path)
#'
#' # label
#' abd <- label_raster_stack(abd)
#' names(abd)
#'
#' }
label_raster_stack <- function(x) {
  stopifnot(inherits(x, "Raster"))

  # check length
  if((raster::nlayers(x) != 52)) {
    stop(paste("The input Raster* object must be full stack or brick of 52",
               "layers as originally provided."))
  }

  srd_date_vec <- seq(from = 0, to = 1, length = 52 + 1)
  srd_date_vec <- (srd_date_vec[1:52] + srd_date_vec[2:(52 + 1)]) / 2
  srd_date_vec <- round(srd_date_vec, digits = 4)

  year_seq <- 2015
  p_time <- strptime(x = paste(round(srd_date_vec * 366), year_seq), "%j %Y")
  date_names <- paste(formatC(p_time$mon + 1, width = 2, format = "d",
                              flag = "0"),
                      formatC(p_time$mday, width = 2, format = "d",
                              flag = "0"),
                      sep = "-")

  names(x) <- paste(2016, date_names, sep = "-")

  return(x)
}


#' Parses the names attached to a Raster from label_raster_stack()
#'
#' The [label_raster_stack()] function labels the dates of the estimate rasters
#' in the format of "XYYYY.MM.DD", because of constraints in the `raster`
#' packge. This function converts that character vector into an ISO compliant
#' Date vector.
#'
#' @param x `Raster` object; full or subset of original eBird Status and
#'   Trends product raster GeoTiff.
#'
#' @return Date vector.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example abundance data
#' sp_path <- download_data("example_data")
#' abd <- load_raster("abundance_umean", sp_path)
#'
#' # label
#' abd <- label_raster_stack(abd)
#'
#' # parse dates
#' parse_raster_dates(abd)
#'
#' }
parse_raster_dates <- function(x) {
  stopifnot(inherits(x, "Raster"))
  if(length(grep("X2016.", names(x)[1])) == 0) {
    stop("Raster names not in correct format, call label_raster_stack() first.")
  }

  as.Date(names(x), format = "X%Y.%m.%d")
}


#' Stixel summary file loader
#'
#' Internal function used by [load_pis()] and [load_pds()] to get the stixel
#' summary information (from summary.txt).
#'
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing stixel summary information about each stixel
#'   centroid.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # stixel summaries
#' summaries <- ebirdst:::load_summary(sp_path)
#' dplyr::glimpse(summaries)
#'
#' }
load_summary <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  stixel_path <- file.path(path, "results", "abund_preds", "unpeeled_folds")
  summary_file <- file.path(stixel_path, "summary.txt")

  if(!file.exists(summary_file)) {
    stop(paste0("The file summary.txt does not exist at ", stixel_path))
  }

  summary_df <- data.table::fread(summary_file,
                                  col.names = ebirdst_summary_names,
                                  #col.names = ebirdst_summary_names_tidy,
                                  stringsAsFactors = FALSE,
                                  showProgress = FALSE)
  summary_df <- summary_df[, c("stixel", "data_type") := NULL]
  summary_df <- summary_df[!is.na(summary_df$lon), ]

  summary_df <- as.data.frame(summary_df)
  if (isTRUE(return_sf)) {
    summary_df <- sf::st_as_sf(summary_df, coords = c("lon", "lat"), crs = 4326)
  }
  return(summary_df)
}


#' Load predictor importances for single species eBird Status and Trends products
#'
#' Loads the predictor importance data (from pi.txt), joins with stixel summary
#' data, sets column names, and cleans up the data.frame.
#'
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing predictor importance values for each stixel, as
#'   well as stixel summary information.
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # load predictor importance
#' pis <- load_pis(sp_path)
#'
#' # plot the top 15 predictor importances
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pis(pis, ext = e, n_top_pred = 15, by_cover_class = TRUE)
#'
#' }
load_pis <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  stixel_path <- file.path(path, "results", "abund_preds", "unpeeled_folds")
  pi_file <- file.path(stixel_path, "pi.txt")

  if(!file.exists(pi_file)) {
    stop(paste("The file pi.txt does not exist at:", stixel_path))
  }

  # pi file
  pi_df <- data.table::fread(pi_file,
                             col.names = ebirdst_pi_names,
                             #col.names = ebirdst_pi_names_tidy,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  pi_df <- pi_df[, c("stixel", "data_type") := NULL]

  # summary file
  summary_df <- load_summary(path)
  summary_df <- summary_df[, 1:12]

  # merge pis with summary
  pi_summary <- dplyr::inner_join(summary_df, pi_df, by = "stixel_id")

  pi_summary <- as.data.frame(pi_summary)
  if (isTRUE(return_sf)) {
    pi_summary <- sf::st_as_sf(pi_summary, coords = c("lon", "lat"), crs = 4326)
  }
  return(pi_summary)
}


#' Load partial dependencies for single species eBird Status and Trends products
#'
#' Loads the partial dependency data (from pd.txt), joins with stixel summary
#' data, sets the names`, and cleans up the data.frame. This is one of the
#' slower functions in the package, due to the size of the pd.txt file (usually
#' multiple GB).
#'
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing partial dependency values for each stixel, as
#'   well as stixel summary information. To make these data more compact,
#'   they're stored in a wide format. Each row corresponds to the partial
#'   dependence relationship of one predictor for a single stixel. There are
#'   50 columns (x1-x50) at which the partial dependence is measured and 50
#'   columns (y1-y50) that give the resulting probability of occurrence on the
#'   logit scale.
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # load partial dependence
#' pds <- load_pds(sp_path)
#'
#' # plot partial dependence for effort hours
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' plot_pds(pds, "effort_hrs", ext = e)
#'
#' }
load_pds <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  stixel_path <- file.path(path, "results", "abund_preds", "unpeeled_folds")
  pd_file <- file.path(stixel_path, "pd.txt")

  if(!file.exists(pd_file)) {
    stop(paste("The file pd.txt does not exist at:", stixel_path))
  }

  # load pd.txt
  pd_df <- data.table::fread(pd_file,
                             col.names = ebirdst_pd_names,
                             #col.names = ebirdst_pi_names_tidy,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  pd_df <- pd_df[, c("stixel", "data_type", "spacer") := NULL]

  # load summary file
  summary_df <- load_summary(path)
  summary_df <- summary_df[, 1:12]

  # merge
  pd_summary <- dplyr::right_join(pd_df, summary_df, by = "stixel_id")
  rm(pd_df, summary_df)

  pd_summary <- as.data.frame(pd_summary)
  if (isTRUE(return_sf)) {
    pd_summary <- sf::st_as_sf(pd_summary, coords = c("lon", "lat"), crs = 4326)
  }
  return(pd_summary)
}


#' Raw test data loader
#'
#' Internal function used by [load_test_data()] to get the full test data for
#' calculating predictive performance metrics. This file contains the observed
#' counts and model predicted relative occurrence and abundance values.
#'
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing test data, including checklist locations,
#'   observed counts, and predicted mean and associated error for relative
#'   occurrence and abundance.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # test data
#' test_data <- ebirdst:::load_test_data_raw(sp_path)
#' dplyr::glimpse(test_data)
#'
#' }
load_test_data_raw <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))

  data_path <- file.path(path, "data")
  td_file <- list.files(data_path,
                        pattern = "erd\\.test\\.data\\.csv",
                        full.names = TRUE)

  if (length(td_file) != 1 || !file.exists(td_file)) {
    stop(paste("Error locating raw test data file at:", data_path))
  }

  # load pd.txt
  td_df <- data.table::fread(td_file,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  td_df <- td_df[, "type" := NULL]

  # fix names
  names(td_df) <- stringr::str_replace_all(stringr::str_to_lower(names(td_df)),
                                           "\\.", "_")
  nm_idx <- match(c("longitude", "latitude"), names(td_df))
  names(td_df)[nm_idx] <- c("lon", "lat")

  td_df <- as.data.frame(td_df)
  if (isTRUE(return_sf)) {
    td_df <- sf::st_as_sf(td_df, coords = c("lon", "lat"), crs = 4326)
  }
  return(td_df)
}


#' Test data loader
#'
#' Internal function used by [compute_ppms()] to get the test data for
#' calculating predictive performance metrics. This file contains the observed
#' counts and model predicted relative occurrence and abundance values.
#'
#' @param path character; full path to the directory containing single species
#'   eBird Status and Trends products.
#' @param return_sf logical; whether to return an [sf] object of spatial points
#'   rather then the default data frame.
#'
#' @return data.frame containing test data, including checklist locations,
#'   observed counts, and predicted mean and associated error for relative
#'   occurrence and abundance.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download example data
#' sp_path <- download_data("example_data")
#'
#' # test data
#' test_data <- load_test_data(sp_path)
#' dplyr::glimpse(test_data)
#'
#' }
load_test_data <- function(path, return_sf = FALSE) {
  stopifnot(dir.exists(path))
  stopifnot(is.logical(return_sf), length(return_sf) == 1)

  stixel_path <- file.path(path, "results", "abund_preds", "unpeeled_folds")
  td_file <- file.path(stixel_path, "test.pred.ave.txt")

  if(!file.exists(td_file)) {
    stop(paste("The file test.pred.ave.txt does not exist at:", stixel_path))
  }

  # load pd.txt
  td_df <- data.table::fread(td_file,
                             col.names= ebirdst_td_names,
                             stringsAsFactors = FALSE,
                             showProgress = FALSE)
  td_df <- td_df[, "data_type" := NULL]


  td_df <- as.data.frame(td_df)
  if (isTRUE(return_sf)) {
    td_df <- sf::st_as_sf(td_df, coords = c("lon", "lat"), crs = 4326)
  }
  return(td_df)
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
