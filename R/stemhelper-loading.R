#' Internal function for transforming st_extent to sinusoidal raster extent
#'
get_sinu_ext <- function(extent) {
  # projection information
  ll <- "+init=epsg:4326"

  sp_ext <- raster::extent(extent$lon.min,
                           extent$lon.max,
                           extent$lat.min,
                           extent$lat.max)
  extllr <- raster::raster(ext = sp_ext)
  extllr[is.na(extllr)] <- 0
  raster::crs(extllr) <- ll

  extsinur <- raster::projectRaster(extllr,
                                    crs = sp::proj4string(
                                      stemhelper::template_raster))

  return(raster::extent(extsinur))
}


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
stack_stem <- function(path,
                       variable,
                       ext = NA,
                       use_analysis_extent = TRUE) {

  poss_var <- c("abundance_ensemble_support", "abundance_lower",
                "abundance_upper", "abundance_umean", "occurrence_umean")
  if(!(variable %in% poss_var)) {
    stop(paste("Selected variable is not one of the following: ",
               paste(poss_var, collapse = ", "), ".", sep = ""))
  }

  fp <- paste(path, "/results/tifs/presentation/", variable, sep = "")

  load_extent <- NULL
  load_extend <- FALSE
  e <- load_config(path)

  if(all(is.na(ext))) {
    if(use_analysis_extent == TRUE) {
      # load with extent
      if(!is.null(e$SPATIAL_EXTENT_LIST)) {
        load_extent <- get_sinu_ext(e$SPATIAL_EXTENT_LIST)
      } else {
        load_extend <- TRUE
      }
    } else {
      if(e$SRD_AGG_FACT == 1) {
        # load without extend
        load_extend <- FALSE
      } else {
        # load with extend
        load_extend <- TRUE
      }
    }
  } else {
    # load with extent
    load_extent <- get_sinu_ext(ext)
  }

  # define function to load and extend each file in path
  load_and_extend <- function(x, ext, use_extend) {
    if(tools::file_ext(x) == "tif") {
      if(!is.null(ext)) {
        r <- raster::crop(
          raster::extend(
            raster::raster(paste(fp, "/", x, sep = "")), ext), ext)
      } else {
        if(use_extend == TRUE) {
          r <- raster::extend(raster::raster(paste(fp, "/", x, sep = "")),
                              stemhelper::template_raster)
        } else {
          r <- raster::raster(paste(fp, "/", x, sep = ""))
        }
      }

      return(r)
    }
  }

  # check to see if path contains more than 1 geotiff file
  if( sum(tools::file_ext(list.files(fp)) == "tif", na.rm = TRUE) < 2 ) {
    stop("Directory does not contain at least 2 .tif files.")
  }

  all_lays <- lapply(X = list.files(fp),
                     FUN = load_and_extend,
                     ext = load_extent,
                     use_extend = load_extend)

  all_lays <- all_lays[!sapply(all_lays, is.null)]

  st <- raster::stack(all_lays)
  rm(all_lays)

  return(st)
}

#' Config file loader
#'
#' Internal function used by load_summary(), load_pis(), and load_pds() to get
#' configuration variables from species run information.
load_config <- function(path) {
  e <- new.env()

  config_file_path <- list.files(paste(path, "/data", sep = ""),
                                 pattern = "*_config*")
  config_file <- paste(path, "/data/", config_file_path, sep = "")

  if(!file.exists(config_file)) {
    stop("*_config.RData file does not exist in the /data directory.")
  }

  load(config_file, envir = e)
  rm(config_file)

  return(e)
}

#' Summary file loader
#'
#' Internal function used by load_pis() and load_pds() to get the stixel
#' summary information
load_summary <- function(path) {
  e <- load_config(path)

  # define stixel summary fields
  train_covariate_means_names <- paste("train.cov.mean",
                                       e$PREDICTOR_LIST,
                                       sep = "_")
  srd_covariate_means_names <- paste("srd.cov.mean",
                                     e$PREDICTOR_LIST,
                                     sep = "_")
  summary_vec_name_vec <- c("srd.n",
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
                            "binary_Kappa",
                            "binary_AUC",
                            "binary.deviance_model",
                            "binary.deviance_mean",
                            "binary.deviance_explained",
                            "pois.deviance_model",
                            "pois.deviance_mean",
                            "posi.deviance_explained",
                            "total_EFFORT_HRS",
                            "total_EFFORT_DISTANCE_KM",
                            "total_NUMBER_OBSERVERS",
                            "train_elevation_mean",
                            train_covariate_means_names, #k-covariate values
                            "train_covariate_entropy",
                            "srd_elevation_mean",
                            srd_covariate_means_names, #k-covariate values
                            "srd_covariate_entropy")

  stixel_path <- "/results/abund_preds/unpeeled_folds/"
  summary_file <- paste(path, stixel_path, "summary.txt", sep = "")

  if(!file.exists(summary_file)) {
    stop(paste("The file summary.txt does not exist at ",
               path, stixel_path, sep = ""))
  }

  summary_vec <- data.table::fread(summary_file)
  names(summary_vec)[3] <- "stixel.id"
  names(summary_vec)[4:ncol(summary_vec)] <- summary_vec_name_vec

  summary_nona <- summary_vec[!is.na(summary_vec$centroid.lon), ]
  rm(summary_vec, summary_file)

  return(summary_nona)
}

#' Load Predictor Importance file for a single species
#'
#' @import data.table
#' @export
load_pis <- function(path) {
  # load config vars
  e <- load_config(path)

  # load pi.txt and set column names
  stixel_path <- "/results/abund_preds/unpeeled_folds/"
  pi_file <- paste(path, stixel_path, "pi.txt", sep = "")

  if(!file.exists(pi_file)) {
    stop(paste("The file pi.txt does not exist at",
               path, stixel_path, sep = ""))
  }

  pi_vec <- data.table::fread(pi_file)
  names(pi_vec)[4:ncol(pi_vec)] <- e$PI_VARS
  names(pi_vec)[3] <- "stixel.id"

  # get summary file
  summary_file <- load_summary(path)

  # merge pis with summary
  pi_summary <- merge(pi_vec, summary_file, by = c("stixel.id"))
  rm(pi_vec, summary_file)

  # return subset
  pi_summary[, c("V1.x", "V2.x", "V1.y", "V2.y") := NULL]
  pi_summary <- pi_summary[, 1:(length(e$PI_VARS) + 12)]

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
  pd_file <- paste(path, stixel_path, "/pd.txt", sep = "")

  if(!file.exists(pd_file)) {
    stop(paste("The file pd.txt does not exist at ",
               path, stixel_path, sep = ""))
  }

  # load pd.txt
  pd_vec <- data.table::fread(pd_file, showProgress = FALSE)
  names(pd_vec)[3] <- "stixel.id"

  # load summary file
  summary_file <- load_summary(path)

  # merge
  pd_summary <- merge(pd_vec, summary_file, by = c("stixel.id"), all.y = TRUE)
  rm(pd_vec, summary_file)

  # return a subset of the fields
  pd_summary[ ,c("V1.x", "V2.x", "V1.y", "V2.y") := NULL]
  pd_summary <- pd_summary[, 1:(length(e$PD_VARS) + 28)]

  return(as.data.frame(pd_summary))
}
