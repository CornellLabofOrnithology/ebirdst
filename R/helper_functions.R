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
#'
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

#' Calculates spatial extent of a stack for plotting
#'
#' After creating a stack, there are lots of NA values and plots of the
#' individual raster layers are at the full extent of the template_raster. To
#' show an ideal extent, this function trims away 0 and NA values and checks
#' to make sure it returns a reasonable extent (based on the template_raster)
#' for plotting. The returned extent object can then be used for plotting.
#'
#' @usage \code{calc_full_extent(stack))}
#'
#' @param stack A RasterStack object, ideally of occurrence or abundace.
#'
#' @return raster Extent object
#'
#' @examples
#' tif_path <- "~"
#' raster_stack <- stack_stem(tif_path)
#' plot_extent <- calc_full_extent(raster_stack)
#' raster::plot(raster_stack[[1]], ext=plot_extent)
calc_full_extent <- function(stack) {

  # aggregate stack for speed, otherwise everything else takes too long
  stack <- raster::aggregate(stack, fact=3)

  # convert 0s to NAs, otherwise trimming is slow and the extent is too broad
  stack[stack == 0] <- NA

  # trim away 0s to get closest extent to positive values
  stack <- raster::trim(stack, values = NA)

  # save extent
  map_extent <- raster::extent(stack)

  # sometimes extent calculations get weird and you'll get a very broad
  # extent that goes further than you want
  # this section is intended to correct that by setting the mins and maxes
  # to that of the template_raster
  # however, if the stack comes in a different projection, we need to correct
  # for that

  # if stack projection is different from template_raster,
  # project template_raster
  if(raster::projection(template_raster) != raster::projection(stack)) {
    this_template_raster <- projectRaster(template_raster,
                                          raster::projection(stack))
  } else {
    this_template_raster <- template_raster
  }

  # create object of this template_raster for comparison to extent extracted
  # from the stack above
  template_raster_extent <- raster::extent(this_template_raster)

  # xmin too low
  if(map_extent[1] < template_raster_extent[1]) {
    map_extent[1] <- template_raster_extent[1]
  }

  # xmax too high
  if(map_extent[2] > template_raster_extent[2]) {
    map_extent[2] <- template_raster_extent[2]
  }

  # ymin too low
  if(map_extent[3] < template_raster_extent[3]) {
    map_extent[3] <- template_raster_extent[3]
  }

  # ymax too high
  if(map_extent[4] > template_raster_extent[4]) {
    map_extent[4] <- template_raster_extent[4]
  }

  return(map_extent)
}

#' Calculates bins based on standard deviations of log-transformed data
#'
#' Mapping species abundnace across the full-annual cycle presents a challenge, in that patterns of concentration and dispersion in abundance change throughout the year, making it difficult to define color bins that suit all seasons and accurately reflect the detail of abundance predictions. To address this, we selected a method (described by Maciejewski et al. 2013) that log transforms the entire year of data, constructs bins with the log-transformed data using standard-deviations, and then untransforms the bins.
#'
#' @usage \code{calc_bins(stack)}
#'
#' @param stack A RasterStack object, of abundance results.
#'
#' @return vector containing break points of bins
#'
#' @examples
#' tif_path <- "~"
#' raster_stack <- stack_stem(tif_path)
#' year_bins <- calc_bins(raster_stack)
#'
#' raster::plot(raster_stack[[1]], xaxt='n', yaxt='n', breaks=year_bins)
calc_bins <- function(stack) {
  zrv <- raster::getValues(stack)
  lzwk <- log(zrv[!is.na(zrv)])
  rm(zrv)

  # LOG SD
  mdl <- mean(lzwk)
  sdl <- sd(lzwk)
  log_sd <- c(mdl-(3.00*sdl),mdl-(2.50*sdl),mdl-(2.00*sdl),mdl-(1.75*sdl),
              mdl-(1.50*sdl),mdl-(1.25*sdl),mdl-(1.00*sdl),mdl-(0.75*sdl),
              mdl-(0.50*sdl),mdl-(0.25*sdl),mdl-(0.125*sdl),
              mdl,
              mdl+(0.125*sdl),mdl+(0.25*sdl),mdl+(0.50*sdl),mdl+(0.75*sdl),
              mdl+(1.00*sdl),mdl+(1.25*sdl),mdl+(1.50*sdl),mdl+(1.75*sdl),
              mdl+(2.00*sdl),mdl+(2.50*sdl),mdl+(3.00*sdl))

  if(max(lzwk) > mdl+(3.00*sdl)) {
    log_sd <- append(log_sd, max(lzwk))
  }

  if(max(lzwk) < mdl+(3.00*sdl)) {
    log_sd <- log_sd[1:length(log_sd)-1]
  }

  if(min(lzwk) < mdl-(3.00*sdl)) {
    log_sd <- append(log_sd, min(lzwk), after=0)
  }

  if(min(lzwk) > mdl-(3.00*sdl)) {
    log_sd <- log_sd[2:length(log_sd)]
  }

  if(exp(min(lzwk)) > 0) {
    log_sd <- append(log_sd, -Inf, after=0)
  }

  rm(lzwk)

  bins <- exp(log_sd)

  return(bins)
}

#' Loads PI data
#'
load_pis <- function(path) {
  # load config vars
  e <- new.env()
  config_file <- list.files(path, pattern="*_config*")
  load(paste(path, "/", config_file, sep = ""), envir = e)

  # load pi.txt
  pi_vec <- data.table::fread(paste(path, "/stixels/pi.txt", sep = ""))
  names(pi_vec)[4:ncol(pi_vec)] <- e$PI_VARS
  names(pi_vec)[3] <- "stixel.id"

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
  summary_vec <- data.table::fread(paste(path, "/stixels/summary.txt", sep=""))
  names(summary_vec)[3] <- "stixel.id"
  names(summary_vec)[4:ncol(summary_vec)] <- summary_vec_name_vec

  summary_nona <- summary_vec[!is.na(summary_vec$centroid.lon), ]

  pi_summary <- merge(pi_vec, summary_nona, by = c("stixel.id"))

  # TODO, what else should we select to return?
  pi_summary[,c("V1.x", "V2.x", "V1.y", "V2.y") := NULL]

  return(pi_summary)
}
