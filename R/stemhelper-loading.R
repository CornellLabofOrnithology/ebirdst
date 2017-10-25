#' Projects st_extent lat/lon list to sinusoidal raster Extent
#'
#' Internal function that converts the st_extent list used throughout this
#' package from lat/lon corners to a raster Extent using the same
#' Sinusoidal projection as the `template_raster` data object.
#'
#' @param st_extent list; st_extent list with lat/lon coordinates.
#'
#' @return A raster Extent in Sinusoidal projection.
#'
#' @examples
#' # define st_extent list
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' # convert
#' sinu_e <- get_sinu_ext(ne_extent)
#' sinu_e
get_sinu_ext <- function(st_extent) {
  # projection information
  ll <- "+init=epsg:4326"

  sp_ext <- raster::extent(st_extent$lon.min,
                           st_extent$lon.max,
                           st_extent$lat.min,
                           st_extent$lat.max)
  extllr <- raster::raster(ext = sp_ext)
  extllr[is.na(extllr)] <- 0
  raster::crs(extllr) <- ll

  extsinur <- raster::projectRaster(extllr,
                                    crs = sp::proj4string(
                                      stemhelper::template_raster))

  return(raster::extent(extsinur))
}

#' Spatiotemporal subsetter for STEM result data objects
#'
#' Internal function that takes a data.frame or SpatialPointsDataFrame and
#' a st_extent list and returns a subset of the data object. Currently
#' designed to handle either a 'rectangle' as defined by a lat/lon bounding
#' box or a 'polygon' as defined by a SpatialPolygon* object. The `use_time`
#' parameter allows for subsetting with or without temporal information. The
#' t.min and t.max objects in the `st_extent` list are currently able to wrap
#' time around the year (e.g., t.min = 0.9 and t.max = 0.1 is acceptable).
#'
#' @param data data.frame or SpatialPointsDataFrame; originating from STEM results.
#' @param st_extent list; st_extent list object
#' @param use_time Boolean; indicating whether to use time in subsetting or not.
#'
#' @return Subset of input data as same type.
#'
#' @examples
#' \dontrun{
#'
#' # define species STEM results location and load pis
#' sp_path <- "path to species STEM results"
#' pis <- load_pis(sp_path)
#'
#' # define st_extent list
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' #subset pis with st_extent list
#' st_extent_subset(data = pis, st_extent = ne_extent, use_time = TRUE)
#' }
st_extent_subset <- function(data, st_extent, use_time = TRUE) {
  if(st_extent$type == "rectangle") {
    if(use_time == TRUE) {
      if(st_extent$t.min > st_extent$t.max) {
        subset_data <- data[(data$date > st_extent$t.min |
                              data$date <= st_extent$t.max) &
                              data$lat > st_extent$lat.min &
                              data$lat <= st_extent$lat.max &
                              data$lon > st_extent$lon.min &
                              data$lon <= st_extent$lon.max, ]
      } else {
        subset_data <- data[data$date > st_extent$t.min &
                              data$date <= st_extent$t.max &
                              data$lat > st_extent$lat.min &
                              data$lat <= st_extent$lat.max &
                              data$lon > st_extent$lon.min &
                              data$lon <= st_extent$lon.max, ]
      }
    } else {
      subset_data <- data[data$lat > st_extent$lat.min &
                          data$lat <= st_extent$lat.max &
                          data$lon > st_extent$lon.min &
                          data$lon <= st_extent$lon.max, ]
    }
  } else if(st_extent$type == "polygon") {
    ll_epsg <- "+init=epsg:4326"
    ll <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

    # check Polygons
    if(is.null(st_extent$polygon)) {
      stop("polygon data not present.")
    }

    # cast data to SpatialPoints
    if(is.data.frame(data)) {
      sites <- sp::SpatialPointsDataFrame(coords = data[, c("lon", "lat")],
                                          data = data,
                                          proj4string = sp::CRS(ll_epsg))
    } else {
      sites <- data
    }

    # check prj
    if(!sp::identicalCRS(sites, st_extent$polygon)) {
      plygn <- sp::spTransform(st_extent$polygon,
                               sp::CRS(sp::proj4string(sites)))
    } else {
      plygn <- st_extent$polygon
    }

    # use for subset
    subset_data <- sites[!is.na(over(sites, as(plygn, "SpatialPolygons"))), ]

    if(is.data.frame(data)) {
      subset_data <- subset_data@data
    }

    if(use_time == TRUE) {
      if(st_extent$t.min > st_extent$t.max) {
        subset_data <- subset_data[(subset_data$date > st_extent$t.min |
                                      subset_data$date <= st_extent$t.max), ]
      } else {
        subset_data <- subset_data[subset_data$date > st_extent$t.min &
                                     subset_data$date <= st_extent$t.max, ]
      }
    }
  } else {
    stop(paste("Spatiotemporal extent type not accepted. ",
               "Use either 'rectangle' or 'polygon'.",
               sep = ""))
  }

  return(subset_data)
}

#' Load, extend, and stack all weeks of STEM rasters for a given result type
#'
#' For one of five result types, loads all available weeks of rasters, extends
#' them to the extent of the study (or to a custom extent), and stacks them
#' into a RasterStack.
#'
#' @param path character; Full path to directory containing the STEM results
#' for a single species.
#' @param variable character; One of: 'abundance_ensemble_support',
#' 'abundance_lower', 'abundance_upper', 'abundance_umean',
#' and 'occurrence_umean'.
#' @param st_extent list; Optional, use to limit the spatial Extent that the rasters are loaded into. Must set use_analysis_extent to FALSE.
#' @param use_analysis_extent Boolean; Default is TRUE. If STEM results were run for a custom non-global extent, that extent object is stored in the configuration file. If TRUE, uses that analysis extent for loading the rasters. If FALSE, rasters are loaded to the full extent of the `template_raster` object.
#'
#' @return RasterStack containing all available weeks of result.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' raster_stack <- stack_stem(path = sp_path, variable = "abundance_umean")
#' }
stack_stem <- function(path,
                       variable,
                       st_extent = NA,
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

  if(all(is.na(st_extent))) {
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
    load_extent <- get_sinu_ext(st_extent)
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
#' configuration variables from STEM species run information
#' (from *_config.RData).
#'
#' @param path character; Full path to the directory containing single species
#' STEM results.
#'
#' @returns environment object containing all run parameters.
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' e <- load_config(sp_path)
#' }
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

#' Stixel summary file loader
#'
#' Internal function used by load_pis() and load_pds() to get the stixel
#' summary information (from summary.txt).
#'
#' @param path character; Full path to the directory containing single species
#' STEM results.
#'
#' @returns data.frame containing stixel summary information about each stixel
#' centroid.
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' summaries <- load_summary(sp_path)
#' }
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
                            "lon",
                            "lat",
                            "date",
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

  summary_vec <- data.table::fread(summary_file, showProgress = FALSE)
  names(summary_vec)[3] <- "stixel.id"
  names(summary_vec)[4:ncol(summary_vec)] <- summary_vec_name_vec

  summary_nona <- summary_vec[!is.na(summary_vec$lon), ]
  rm(summary_vec, summary_file)

  return(summary_nona)
}

#' Load predictor importances for single species STEM results
#'
#' Loads the predictor importance data (from pi.txt), joins with stixel summary
#' data, sets the names from `load_config()`, and cleans up the data.frame.
#'
#' @param path character; Full path to the directory containing single species
#' STEM results.
#'
#' @returns data.frame containing predictor importance values for each stixel,
#' as well as stixel summary information.
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' pis <- load_pis(sp_path)
#' }
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

  pi_vec <- data.table::fread(pi_file, showProgress = FALSE)
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
