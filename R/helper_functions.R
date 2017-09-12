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
#' @export
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
#' @export
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

#' Summary file loader
#'
#' Used by load_pis() and load_pds()
load_summary <- function(path) {
  e <- new.env()
  config_file <- list.files(path, pattern="*_config*")
  load(paste(path, "/", config_file, sep = ""), envir = e)

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

  return(summary_nona)
}

#' Loads PI data
#'
#' @import data.table
#' @export
load_pis <- function(path) {
  # load config vars
  e <- new.env()
  config_file <- list.files(path, pattern="*_config*")
  load(paste(path, "/", config_file, sep = ""), envir = e)

  # load pi.txt
  pi_vec <- data.table::fread(paste(path, "/stixels/pi.txt", sep = ""))
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
  e <- new.env()
  config_file <- list.files(path, pattern="*_config*")
  load(paste(path, "/", config_file, sep = ""), envir = e)

  # load pi.txt
  pd_vec <- data.table::fread(paste(path, "/stixels/pd.txt", sep = ""))
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

#' Plot PIs and PDs
#'
#' @export
plot_centroids <- function(pis,
                           pds,
                           st_extent = NA,
                           plot_pis = TRUE,
                           plot_pds = TRUE,
                           ...) {

  if(plot_pds == FALSE & plot_pis == FALSE) {
    stop("Plotting of both PIs and PDs set to FALSE. Nothing to plot!")
  }

  ll <- "+init=epsg:4326"
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

  par(mar=c(0,0,3,0))
  title_text <- ""

  if(plot_pds == TRUE) {
    tpds <- unique(pds[, c("centroid.lon","centroid.lat","centroid.date")])
    tpds_sp <- sp::SpatialPointsDataFrame(tpds[,c("centroid.lon",
                                                  "centroid.lat")],
                                          tpds,
                                          proj4string = sp::CRS(ll))
    tpds_ext <- raster::extent(tpds_sp)
    tpds_prj <- sp::spTransform(tpds_sp, sp::CRS(mollweide))
    rm(tpds)

    # this is the wrong way to check
    if(!is.na(st_extent)) {
      tpds_sub <- tpds_sp[tpds_sp$centroid.date > st_extent$t.min &
                          tpds_sp$centroid.date <= st_extent$t.max &
                          tpds_sp$centroid.lat > st_extent$y.min &
                          tpds_sp$centroid.lat <= st_extent$y.max &
                          tpds_sp$centroid.lon > st_extent$x.min &
                          tpds_sp$centroid.lon <= st_extent$x.max, ]

      tpds_region <- sp::spTransform(tpds_sub, sp::CRS(mollweide))
    }

    rm(tpds_sp)

    # start plot with all possible PDs
    raster::plot(tpds_prj,
                 ext = raster::extent(tpds_prj),
                 col = "darkblue",
                 cex = 0.5,
                 pch = 16)

    # plot PDs in st_extent
    if(!is.na(st_extent)) {
      raster::plot(tpds_region,
                   ext = raster::extent(tpds_prj),
                   col = "blue",
                   cex = 0.5,
                   pch = 16,
                   add = TRUE)
    }

    title_text <- paste(title_text,
                        "Available PDs: ", nrow(tpds_prj), "\n",
                        "Selected PDs: ", nrow(tpds_region), "\n",
                        sep = "")
  }

  if(plot_pis == TRUE) {
    tpis <- unique(pis[, c("centroid.lon","centroid.lat","centroid.date")])
    tpis_sp <- sp::SpatialPointsDataFrame(tpis[,c("centroid.lon",
                                                  "centroid.lat")],
                                          tpis,
                                          proj4string = sp::CRS(ll))
    tpis_ext <- raster::extent(tpis_sp)
    tpis_prj <- sp::spTransform(tpis_sp, mollweide)
    rm(tpis)

    # this is the wrong way to check
    if(!is.na(st_extent)) {
      tpis_sub <- tpis_sp[tpis_sp$centroid.date > st_extent$t.min &
                          tpis_sp$centroid.date <= st_extent$t.max &
                          tpis_sp$centroid.lat > st_extent$y.min &
                          tpis_sp$centroid.lat <= st_extent$y.max &
                          tpis_sp$centroid.lon > st_extent$x.min &
                          tpis_sp$centroid.lon <= st_extent$x.max, ]

      tpis_region <- sp::spTransform(tpis_sub, sp::CRS(mollweide))
    }
    rm(tpis_sp)

    if(plot_pds == TRUE) {
      # use extent from all PDs
      e <- as(tpds_ext, "SpatialPolygons")
    } else {
      # use extent from all PIs
      e <- as(tpis_ext, "SpatialPolygons")
    }

    sp::proj4string(e) <- sp::CRS(sp::proj4string(rnaturalearthdata::countries50))
    c50 <- rnaturalearthdata::countries50
    s50 <- rnaturalearthdata::states50
    wh <- c50[!is.na(sp::over(c50, e)),]
    wh_states <- s50[!is.na(sp::over(s50, e)),]
    wh_moll <- sp::spTransform(wh, sp::CRS(mollweide))
    wh_states_moll <- sp::spTransform(wh_states, sp::CRS(mollweide))

    if(plot_pds == TRUE) {
      # start plot with all possible PDs
      raster::plot(tpis_prj,
                   ext = raster::extent(tpds_prj),
                   col = "darkred",
                   cex = 0.5,
                   pch = 16,
                   add = TRUE)

      # plot PDs in st_extent
      if(!is.na(st_extent)) {
        raster::plot(tpis_region,
                     ext = raster::extent(tpds_prj),
                     col = "red",
                     cex = 0.5,
                     pch = 16,
                     add = TRUE)
      }

      raster::plot(wh_moll, ext = raster::extent(tpds_prj), add = TRUE)
      raster::plot(wh_states_moll,
                   ext = raster::extent(tpds_prj),
                   lwd = 0.5,
                   add = TRUE)
    } else {
      # start plot with all possible PDs
      raster::plot(tpis_prj,
                   ext = raster::extent(tpis_sp_prj),
                   col = "darkred",
                   cex = 0.5,
                   pch = 16)

      # plot PDs in st_extent
      if(!is.na(st_extent)) {
        raster::plot(tpis_region,
                     ext = raster::extent(tpis_sp_prj),
                     add = TRUE,
                     col = "red",
                     cex = 0.5,
                     pch = 2,
                     add = TRUE)
      }
      raster::plot(wh_moll, ext = raster::extent(tpis_sp_prj), add = TRUE)
      raster::plot(wh_states_moll,
                   ext = raster::extent(tpis_sp_prj),
                   lwd = 0.5,
                   add = TRUE)
    }

    title_text <- paste(title_text,
                        "Available PIs: ", nrow(tpis_prj), "\n",
                        "Selected PIs: ", nrow(tpis_region), "\n",
                        sep = "")
  }

  title(main = title_text, cex.main = 0.75, line=-1)
}

#' Plot extent of estimation calculated from subset of centroids
#'
#' @export
#' @import sp
calc_effective_extent <- function(st_extent,
                                  pis = NA,
                                  pds = NA) {

  if(is.na(pis) & is.na(pds)) {
    stop("Both PIs and PDs are NA. Nothing to calculate.")
  }

  if(!is.na(pis) & !is.na(pds)) {
    stop("Unable to calculate for both PIs and PDs, supply one or the other.")
  }


  # select the centroid locs
  ll <- "+init=epsg:4326"

  if(!is.na(pis)) {
    stixels <- pis
  } else {
    stixels <- pds
  }

  tpis <- unique(stixels[, c("centroid.lon", "centroid.lat", "centroid.date",
                         "stixel_width", "stixel_height")])
  tpis_sp <- sp::SpatialPointsDataFrame(tpis[,c("centroid.lon",
                                                "centroid.lat")],
                                        tpis,
                                        proj4string = sp::CRS(ll))

  tpis_sp_prj <- sp::spTransform(tpis_sp,
                                 sp::CRS(sp::proj4string(template_raster)))

  sp_ext <- raster::extent(tpis_sp_prj)
  rm(tpis, tpis_sp_prj)

  tpis_sub <- tpis_sp[tpis_sp$centroid.date > st_extent$t.min &
                        tpis_sp$centroid.date <= st_extent$t.max &
                        tpis_sp$centroid.lat > st_extent$y.min &
                        tpis_sp$centroid.lat <= st_extent$y.max &
                        tpis_sp$centroid.lon > st_extent$x.min &
                        tpis_sp$centroid.lon <= st_extent$x.max, ]

  # build stixels as polygons

  # create corners
  xPlus <- tpis_sub$centroid.lon + (tpis_sub$stixel_width/2)
  yPlus <- tpis_sub$centroid.lat + (tpis_sub$stixel_height/2)
  xMinus <- tpis_sub$centroid.lon - (tpis_sub$stixel_width/2)
  yMinus <- tpis_sub$centroid.lat - (tpis_sub$stixel_height/2)

  ID <- row.names(tpis_sub)

  square <- cbind(xMinus, yPlus, xPlus, yPlus, xPlus,
                  yMinus, xMinus, yMinus, xMinus, yPlus)

  polys <- sp::SpatialPolygons(mapply(function(poly, id) {
    xy <- matrix(poly, ncol=2, byrow=TRUE)
    sp::Polygons(list(sp::Polygon(xy)), ID=id)
  },
  split(square, row(square)),
  ID),
  proj4string=sp::CRS("+init=epsg:4326"))

  tdspolydf <- sp::SpatialPolygonsDataFrame(polys, tpis_sub@data)

  # assign value
  tdspolydf$weight <- 1

  # project to template raster
  tdspolydf_prj <- sp::spTransform(tdspolydf,
                                   sp::CRS(sp::proj4string(template_raster)))

  # summarize...not sure how to do this step
  tpis_r <- raster::rasterize(tdspolydf_prj,
                              template_raster,
                              field="weight",
                              fun=sum)

  tpis_per <- tpis_r/nrow(tpis_sub)
  tpis_per[tpis_per < 0.50] <- NA

  # plot
  countries50_prj <- sp::spTransform(rnaturalearthdata::countries50,
                                      sp::CRS(sp::proj4string(template_raster)))
  states50_prj <- sp::spTransform(rnaturalearthdata::states50,
                                 sp::CRS(sp::proj4string(template_raster)))
  e <- as(sp_ext, "SpatialPolygons")
  sp::proj4string(e) <- sp::CRS(sp::proj4string(template_raster))
  wh <- countries50_prj[!is.na(sp::over(countries50_prj, e)),]
  wh_states <- states50_prj[!is.na(sp::over(states50_prj, e)),]
  mollweide <- sp::CRS("+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84")
  wh_moll <- sp::spTransform(wh, mollweide)
  wh_states_moll <- sp::spTransform(wh_states, mollweide)

  tpis_per_prj <- raster::projectRaster(raster::mask(tpis_per, template_raster),
                                        crs=mollweide)

  tdspolydf_moll <- sp::spTransform(tdspolydf_prj, mollweide)

  # project the selected points to mollweide
  tpis_sub_moll <- sp::spTransform(tpis_sub, mollweide)


  par(mar=c(0,0,0,2))
  raster::plot(tpis_per_prj,
               xaxt='n',
               yaxt='n',
               bty='n',
               ext=raster::extent(tdspolydf_moll),
               col=viridis::viridis(100),
               maxpixels=raster::ncell(tpis_per_prj),
               legend=TRUE)
  sp::plot(wh_moll, add=TRUE, border='gray')
  sp::plot(wh_states_moll, add=TRUE, border='gray')
  sp::plot(tpis_sub_moll, add=TRUE, pch=19)

  return(tpis_per)
}
