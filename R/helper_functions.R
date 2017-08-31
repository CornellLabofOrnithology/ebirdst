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

  rm(lzwk)

  bins <- exp(log_sd)

  return(bins)
}
