#' Calculates spatial extent of Raster* object for plotting
#'
#' After creating a stack, there are lots of NA values and plots of the
#' individual raster layers are at the full extent of the template_raster. To
#' show an ideal extent, this function trims away 0 and NA values and checks
#' to make sure it returns a reasonable extent (based on the template_raster)
#' for plotting. The returned extent object can then be used for plotting.
#'
#' @usage \code{calc_full_extent(x))}
#'
#' @param stack Raster* object, ideally of occurrence or abundace.
#'
#' @return raster Extent object
#' @export
#' @examples
#' tif_path <- "~"
#' raster_stack <- stack_stem(tif_path)
#' plot_extent <- calc_full_extent(raster_stack)
#' raster::plot(raster_stack[[1]], ext = plot_extent)
calc_full_extent <- function(x) {

  # aggregate stack for speed, otherwise everything else takes too long
  stack <- raster::aggregate(x, fact = 3)

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
  tr <- stemhelper::template_raster

  if(raster::projection(tr) != raster::projection(stack)) {
    this_template_raster <- projectRaster(tr, crs = sp::proj4string(stack))
  } else {
    this_template_raster <- tr
  }
  rm(tr, stack)

  # create object of this template_raster for comparison to extent extracted
  # from the stack above
  template_raster_extent <- raster::extent(this_template_raster)
  rm(this_template_raster)

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
#' @usage \code{calc_bins(x)}
#'
#' @param x Raster* object
#'
#' @return vector containing break points of bins
#' @export
#' @examples
#' tif_path <- "~"
#' raster_stack <- stack_stem(tif_path)
#' year_bins <- calc_bins(raster_stack)
#'
#' raster::plot(raster_stack[[1]], xaxt = 'n', yaxt = 'n', breaks = year_bins)
calc_bins <- function(x) {

  # get a vector of all the values in the stack
  zrv <- raster::getValues(x)
  zrv[zrv == 0] <- NA

  # log transform the non-NA values
  # there shouldn't be zero values in here based on how the data is stored
  # however, it might be worth a check
  lzwk <- log(zrv[!is.na(zrv)])
  rm(zrv)

  # setup the binning structure
  # calculate metrics
  maxl <- max(lzwk)
  minl <- min(lzwk)
  mdl <- mean(lzwk)
  sdl <- sd(lzwk)
  rm(lzwk)

  # build a vector of bins
  log_sd <- c(mdl - (3.00 * sdl),
              mdl - (2.50 * sdl),
              mdl - (2.00 * sdl),
              mdl - (1.75 * sdl),
              mdl - (1.50 * sdl),
              mdl - (1.25 * sdl),
              mdl - (1.00 * sdl),
              mdl - (0.75 * sdl),
              mdl - (0.50 * sdl),
              mdl - (0.25 * sdl),
              mdl - (0.125 * sdl),
              mdl,
              mdl + (0.125 * sdl),
              mdl + (0.25 * sdl),
              mdl + (0.50 * sdl),
              mdl + (0.75 * sdl),
              mdl + (1.00 * sdl),
              mdl + (1.25 * sdl),
              mdl + (1.50 * sdl),
              mdl + (1.75 * sdl),
              mdl + (2.00 * sdl),
              mdl + (2.50 * sdl),
              mdl + (3.00 * sdl))

  # lots of checks for values outside of the upper and lower bounds

  # remove +3 SD break if it is greater than max
  if(maxl < mdl + (3.00 * sdl)) {
    log_sd <- log_sd[1:length(log_sd)-1]
  }

  # add max if the max is greater than +3 SD break
  if(maxl > mdl + (3.00 * sdl) | maxl > log_sd[length(log_sd)]) {
    log_sd <- append(log_sd, maxl)
  }

  # remove the -3 SD break if it is less than the min
  if(minl > mdl - (3.00 * sdl)) {
    log_sd <- log_sd[2:length(log_sd)]
  }

  # add min if the min is less than -3 SD break
  if(minl < mdl - (3.00 * sdl) | minl < log_sd[1]) {
    log_sd <- append(log_sd, minl, after = 0)
  }

  # if the untransformed min is greater than 0, add a zero break
  if(exp(minl) > 0) {
    log_sd <- append(log_sd, -Inf, after = 0)
  }

  # untransform
  bins <- exp(log_sd)
  rm(log_sd)

  return(bins)
}

#' Combine zero layers and abundance to create a single layer
#'
#' @export
combine_layers <- function(stack, path, week) {
  e <- load_config(path)

  # possible checks
  # week is between 1 and 52

  # load abundance
  abund_week <- stack[[week]]

  # load positive ensemble support
  pos_es_stack <- stack_stem(path, variable = "abundance_ensemble_support")
  pos_es_week <- pos_es_stack[[week]]

  # load zero ensemble support
  zero_es_week <- zero_es_stack[[week]]

  if(e$SRD_AGG_FACT == 1) {
    zero_es_week <- raster::resample(zero_es_week, pos_es_week)
    zero_es_mask <- raster::mask(zero_es_week, pos_es_week)
  } else {
    zero_es_mask <- raster::mask(zero_es_week, stemhelper::template_raster)
  }

  # set zeroes
  pos_es_week[pos_es_week] <- 0
  zero_es_mask <- zero_es_mask >= 99
  zero_es_mask[zero_es_mask == 0] <- NA
  zero_es_mask[zero_es_mask == 1] <- 0

  # stack and max
  week_stack <- raster::stack(pos_es_week, zero_es_mask, abund_week)
  week_max <- raster::calc(week_stack, max, na.rm = TRUE)

  return(week_max)
}

#' Map PI and PD centroid locations for species
#'
#' @export
map_centroids <- function(pis,
                          pds,
                          st_extent = NA,
                          plot_pis = TRUE,
                          plot_pds = TRUE,
                          ...) {

  if(plot_pds == FALSE & plot_pis == FALSE) {
    stop("Plotting of both PIs and PDs set to FALSE. Nothing to plot!")
  }

  # projection information
  ll <- "+init=epsg:4326"
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

  # initialize graphical parameters
  par(mfrow = c(1, 1), mar=c(0,0,0,0), bg = "black")

  # Plotting PDs
  if(plot_pds == TRUE) {
    tpds <- unique(pds[, c("centroid.lon", "centroid.lat", "centroid.date")])
    tpds_sp <- sp::SpatialPointsDataFrame(tpds[, c("centroid.lon",
                                                   "centroid.lat")],
                                          tpds,
                                          proj4string = sp::CRS(ll))
    tpds_prj <- sp::spTransform(tpds_sp, sp::CRS(mollweide))
    rm(tpds)

    # start plot with all possible PDs
    raster::plot(tpds_prj,
                 ext = raster::extent(tpds_prj),
                 col = "#1b9377",
                 cex = 0.01,
                 pch = 16)

    sp::plot(ned_wh_co_moll, col = "#5a5a5a", add = TRUE)

    raster::plot(tpds_prj,
                 ext = raster::extent(tpds_prj),
                 col = "#1b9377",
                 cex = 0.4,
                 pch = 16,
                 add = TRUE)

    if(!all(is.na(st_extent))) {
      tpds_sub <- tpds_sp[tpds_sp$centroid.date > st_extent$t.min &
                          tpds_sp$centroid.date <= st_extent$t.max &
                          tpds_sp$centroid.lat > st_extent$lat.min &
                          tpds_sp$centroid.lat <= st_extent$lat.max &
                          tpds_sp$centroid.lon > st_extent$lon.min &
                          tpds_sp$centroid.lon <= st_extent$lon.max, ]

      tpds_region <- sp::spTransform(tpds_sub, sp::CRS(mollweide))
      rm(tpds_sub)

      # plot PDs in st_extent
      raster::plot(tpds_region,
                   ext = raster::extent(tpds_prj),
                   col = "#b3e2cd",
                   cex = 0.4,
                   pch = 16,
                   add = TRUE)
    }
    rm(tpds_sp)

    # xmin, xmax, ymin, ymax
    usr <- par("usr")
    xwidth <- usr[2] - usr[1]
    yheight <- usr[4] - usr[3]

    text(x = usr[1] + xwidth/8,
         y = usr[3] + yheight/7,
         paste("Available PDs: ", nrow(tpds_prj), sep = ""),
         cex = 1,
         col = "#1b9377")

    text(x = usr[1] + xwidth/8,
         y = usr[3] + yheight/9,
         paste("Selected PDs: ", nrow(tpds_region), sep = ""),
         cex = 1,
         col = "#b3e2cd")

    wh_extent <- raster::extent(tpds_prj)
    rm(tpds_prj, tpds_region)
  }

  # Plotting PIs
  if(plot_pis == TRUE) {
    tpis <- unique(pis[, c("centroid.lon", "centroid.lat", "centroid.date")])
    tpis_sp <- sp::SpatialPointsDataFrame(tpis[, c("centroid.lon",
                                                   "centroid.lat")],
                                          tpis,
                                          proj4string = sp::CRS(ll))
    tpis_prj <- sp::spTransform(tpis_sp, mollweide)
    rm(tpis)

    if(plot_pds == FALSE) {
      wh_extent <- raster::extent(tpis_prj)

      # start plot with all possible PDs
      raster::plot(tpis_prj,
                   ext = wh_extent,
                   col = "#1b9377",
                   cex = 0.01,
                   pch = 16)

      sp::plot(ned_wh_co_moll, col = "#5a5a5a", add = TRUE)
    }

    # start plot with all possible PIs
    raster::plot(tpis_prj,
                 ext = wh_extent,
                 col = "#d95f02",
                 cex = 0.4,
                 pch = 16,
                 add = TRUE)

    if(!all(is.na(st_extent))) {
      tpis_sub <- tpis_sp[tpis_sp$centroid.date > st_extent$t.min &
                          tpis_sp$centroid.date <= st_extent$t.max &
                          tpis_sp$centroid.lat > st_extent$lat.min &
                          tpis_sp$centroid.lat <= st_extent$lat.max &
                          tpis_sp$centroid.lon > st_extent$lon.min &
                          tpis_sp$centroid.lon <= st_extent$lon.max, ]

      tpis_region <- sp::spTransform(tpis_sub, sp::CRS(mollweide))

      # plot PIs in st_extent
      raster::plot(tpis_region,
                   ext = wh_extent,
                   col = "#fdcdac",
                   cex = 0.4,
                   pch = 16,
                   add = TRUE)
    }
    rm(tpis_sp)

    # xmin, xmax, ymin, ymax
    usr <- par("usr")
    xwidth <- usr[2] - usr[1]
    yheight <- usr[4] - usr[3]

    text(x = usr[1] + xwidth/8,
         y = usr[3] + yheight/12,
         paste("Available PIs: ", nrow(tpis_prj), sep = ""),
         cex = 1,
         col = "#d95f02")

    text(x = usr[1] + xwidth/8,
         y = usr[3] + yheight/17,
         paste("Selected PIs: ", nrow(tpis_region), sep = ""),
         cex = 1,
         col = "#fdcdac")

    rm(tpis_prj, tpis_region)
  }

  # plot reference data
  raster::plot(ned_wh_co_moll,
               ext = wh_extent,
               lwd = 0.25,
               border = 'black',
               add = TRUE)

  raster::plot(ned_wh_st_moll,
               ext = wh_extent,
               lwd = 0.25,
               border = 'black',
               add = TRUE)
}

#' Map extent of estimation calculated from subset of centroids
#'
#' @export
#' @import sp
calc_effective_extent <- function(st_extent,
                                  pis = NA,
                                  pds = NA) {

  if(is.null(nrow(pis)) & is.null(nrow(pds))) {
    stop("Both PIs and PDs are NA. Nothing to calculate.")
  }

  if(!is.null(nrow(pis)) & !is.null(nrow(pds))) {
    stop("Unable to calculate for both PIs and PDs, supply one or the other.")
  }

  # projection info
  ll <- "+init=epsg:4326"
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

  # set object based on whether using PIs or PDs
  if(!is.null(nrow(pis))) {
    stixels <- pis
  } else {
    stixels <- pds
  }
  rm(pis, pds)

  # subset, create spatial data, project
  tpis <- unique(stixels[, c("centroid.lon", "centroid.lat", "centroid.date",
                             "stixel_width", "stixel_height")])
  tpis_sp <- sp::SpatialPointsDataFrame(tpis[,c("centroid.lon",
                                                "centroid.lat")],
                                        tpis,
                                        proj4string = sp::CRS(ll))
  rm(tpis)

  tpis_sub <- tpis_sp[tpis_sp$centroid.date > st_extent$t.min &
                      tpis_sp$centroid.date <= st_extent$t.max &
                      tpis_sp$centroid.lat > st_extent$lat.min &
                      tpis_sp$centroid.lat <= st_extent$lat.max &
                      tpis_sp$centroid.lon > st_extent$lon.min &
                      tpis_sp$centroid.lon <= st_extent$lon.max, ]
  rm(tpis_sp)

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
    xy <- matrix(poly, ncol = 2, byrow = TRUE)
    sp::Polygons(list(sp::Polygon(xy)), ID = id)
  }, split(square, row(square)), ID), proj4string = sp::CRS("+init=epsg:4326"))

  tdspolydf <- sp::SpatialPolygonsDataFrame(polys, tpis_sub@data)
  rm(xPlus, yPlus, xMinus, yMinus, ID, square, polys)

  # assign value
  tdspolydf$weight <- 1

  # project to template raster
  tdspolydf_prj <- sp::spTransform(tdspolydf,
                                   sp::CRS(
                                     sp::proj4string(
                                       stemhelper::template_raster)))
  rm(tdspolydf)

  # summarize...not sure how to do this step
  tpis_r <- raster::rasterize(tdspolydf_prj,
                              stemhelper::template_raster,
                              field="weight",
                              fun=sum)

  tpis_per <- tpis_r/nrow(tpis_sub)
  rm(tpis_r)
  tpis_per[tpis_per < 0.50] <- NA

  # plot
  tpis_per_prj <- raster::projectRaster(
    raster::mask(tpis_per, stemhelper::template_raster), crs = mollweide)

  tdspolydf_moll <- sp::spTransform(tdspolydf_prj, mollweide)

  # project the selected points to mollweide
  tpis_sub_moll <- sp::spTransform(tpis_sub, mollweide)

  par(mar=c(0,0,0,2))
  raster::plot(tpis_per_prj,
               xaxt = 'n',
               yaxt = 'n',
               bty = 'n',
               ext = raster::extent(tdspolydf_moll),
               col = viridis::viridis(100),
               maxpixels = raster::ncell(tpis_per_prj),
               legend = TRUE)
  sp::plot(ned_wh_co_moll, add = TRUE, border = 'gray')
  sp::plot(ned_wh_st_moll, add = TRUE, border = 'gray')
  sp::plot(tpis_sub_moll, add = TRUE, pch = 16, cex = 1 * par()$cex)

  return(tpis_per)
}
