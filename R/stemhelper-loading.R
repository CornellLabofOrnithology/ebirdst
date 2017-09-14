#' Load, extend, and stack all STEM .tif rasters in a directory
#'
#' Takes all of the .tif rasters in a directory at provided path, loads them,
#' extends them to the extent of the study, and stacks them into a RasterStack.
#' In practice, this will often be all 52 weeks of a single variable
#' (e.g., abundance_mean), but the files could be rearranged and stacked as well
#' (i.e., a single week of a single variable across multiple species).
#'
#' @usage \code{stack_stem(path)}
#'
#' @param path Full path to directory containing more than one STEM .tif raster.
#'
#' @return RasterStack object
#' @export
#' @examples
#' tif_path <- "~"
#' raster_stack <- stack_stem(tif_path)
stack_stem <- function(path) {

  # define function to load and extend each file in path
  load_and_extend <- function(x) {

    if(tools::file_ext(x) == "tif") {
      r <- raster::extend(raster::raster(paste(path, "/", x, sep = "")),
                          template_raster)

      return(r)
    }
  }

  # check to see if path contains more than 1 geotiff file
  if( sum(tools::file_ext(list.files(path)) == "tif", na.rm = TRUE) < 2 ) {
    stop("Directory does not contain at least 2 .tif files.")
  }

  all_lays <- lapply(X = list.files(path), FUN = load_and_extend)

  st <- raster::stack(all_lays)
  rm(all_lays)

  return(st)
}

#' Config file loader
#'
#' Used by load_summary(), load_pis(), and load_pds()
load_config <- function(path) {
  e <- new.env()
  config_file <- list.files(paste(path, "/data", sep = ""), pattern="*_config*")
  load(paste(path, "/data/", config_file, sep = ""), envir = e)

  return(e)
}

#' Summary file loader
#'
#' Used by load_pis() and load_pds()
load_summary <- function(path) {
  e <- load_config(path)

  # load summary
  train_covariate_means_names <- paste("train.cov.mean",
                                       e$PREDICTOR_LIST,
                                       sep = "_")
  srd_covariate_means_names <- paste("srd.cov.mean",
                                     e$PREDICTOR_LIST,
                                     sep = "_")
  summary_vec_name_vec <-c(
    "srd.n",
    "centroid.lon",
    "centroid.lat",
    "centroid.date",
    "stixel_width",
    "stixel_height",
    "stixel_area",
    "train.n",
    "positive.ob_n",
    "stixel_prevalence",
    "mean_non_zero_count",
    # ------------
    "binary_Kappa",
    "binary_AUC",
    # ------------
    "binary.deviance_model",
    "binary.deviance_mean",
    "binary.deviance_explained",
    "pois.deviance_model",
    "pois.deviance_mean",
    "posi.deviance_explained",
    # ------------
    "total_EFFORT_HRS",
    "total_EFFORT_DISTANCE_KM",
    "total_NUMBER_OBSERVERS",
    "train_elevation_mean",
    # ------------
    train_covariate_means_names, #k-covariate values
    "train_covariate_entropy",
    "srd_elevation_mean",
    srd_covariate_means_names, #k-covariate values
    "srd_covariate_entropy" )

  stixel_path <- "/results/abund_preds/unpeeled_folds/"
  summary_vec <- data.table::fread(paste(path, stixel_path,
                                         "summary.txt", sep=""))
  names(summary_vec)[3] <- "stixel.id"
  names(summary_vec)[4:ncol(summary_vec)] <- summary_vec_name_vec

  summary_nona <- summary_vec[!is.na(summary_vec$centroid.lon), ]

  return(summary_nona)
}

#' Loads PI data
#'
#' @import data.table
#' @export
load_pis <- function(path) {
  # load config vars
  e <- load_config(path)

  # load pi.txt
  stixel_path <- "/results/abund_preds/unpeeled_folds/"
  pi_vec <- data.table::fread(paste(path, stixel_path, "pi.txt", sep = ""))
  names(pi_vec)[4:ncol(pi_vec)] <- e$PI_VARS
  names(pi_vec)[3] <- "stixel.id"

  # get summary file
  summary_file <- stemhelper::load_summary(path)

  # merge
  pi_summary <- merge(pi_vec, summary_file, by = c("stixel.id"))
  rm(pi_vec, summary_file)

  # return
  # TODO, what else should we select to return?
  pi_summary[,c("V1.x", "V2.x", "V1.y", "V2.y") := NULL]
  pi_summary <- pi_summary[,1:98]

  return(as.data.frame(pi_summary))
}

#' Loads PD data
#'
#' @import data.table
#' @export
load_pds <- function(path) {
  # load config vars
  e <- load_config(path)

  # load pi.txt
  stixel_path <- "/results/abund_preds/unpeeled_folds/"
  pd_vec <- data.table::fread(paste(path, stixel_path, "/pd.txt", sep = ""))
  names(pd_vec)[3] <- "stixel.id"

  # load summary file
  summary_file <- load_summary(path)

  # merge
  pd_summary <- merge(pd_vec, summary_file, by = c("stixel.id"), all.y = TRUE)
  rm(pd_vec, summary_file)

  # return
  # TODO, what else should we select to return?
  pd_summary[,c("V1.x", "V2.x", "V1.y", "V2.y") := NULL]
  pd_summary <- pd_summary[,1:114]

  return(as.data.frame(pd_summary))
}
