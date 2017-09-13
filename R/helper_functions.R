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

    # plot
    wh <- rnaturalearth::ne_countries(continent = c("North America",
                                                    "South America"),
                                      scale = 50)
    wh_states <- rnaturalearth::ne_states(iso_a2 = unique(wh@data$iso_a2))
    wh_moll <- sp::spTransform(wh, mollweide)
    wh_states_moll <- sp::spTransform(wh_states, mollweide)

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
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

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
  wh <- rnaturalearth::ne_countries(continent = c("North America",
                                                  "South America"),
                                    scale = 50)
  wh_states <- rnaturalearth::ne_states(iso_a2 = unique(wh@data$iso_a2))
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

#' Plotting PIs as barplots
#'
#' @export
#' @import ggplot2 mgcv
plot_pis <- function(pis,
                     st_extent,
                     by_cover_class = FALSE,
                     num_top_preds = 50) {

  # subset for exetnt
  ttt <- pis[pis$centroid.date > st_extent$t.min &
             pis$centroid.date <= st_extent$t.max &
             pis$centroid.lat > st_extent$y.min &
             pis$centroid.lat <= st_extent$y.max &
             pis$centroid.lon > st_extent$x.min &
             pis$centroid.lon <= st_extent$x.max, 2:87]

  if(by_cover_class == TRUE) {

    land.cover.class.codes <- c(1:10,12,13,16)
    lc.tag <- "UMD_FS_C"
    water.cover.class.codes <- c(0,2,3,5,6,7)
    wc.tag <- "MODISWATER_FS_C"
    predictor.names <- names(tpis_sub)

    # ---
    ttt.new <- NULL
    for (iii.pred in 1:length(land.cover.class.codes)){
      # iii.pred <- 1
      new.pred.name <- paste(lc.tag,land.cover.class.codes[iii.pred],"_",sep="")
      pred.nindex <- grep(
        new.pred.name,
        x = predictor.names)
      ttt.new <- cbind(ttt.new,
                       apply( ttt[, pred.nindex], 1, mean, na.rm=T))
      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <-
        paste(lc.tag,land.cover.class.codes[iii.pred],sep="")
    }
    for (iii.pred in 1:length(water.cover.class.codes)){
      # iii.pred <- 1
      new.pred.name <- paste(wc.tag,water.cover.class.codes[iii.pred],"_",sep="")
      pred.nindex <- grep(
        new.pred.name,
        x = predictor.names)
      ttt.new <- cbind(ttt.new,
                       apply( ttt[, pred.nindex], 1, mean, na.rm=T))
      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <-
        paste(wc.tag,water.cover.class.codes[iii.pred],sep="")
    }
    #head(ttt.new)
    ttt <- ttt.new
  }

  # compute median
  pi_median <- apply(ttt, 2, median, na.rm = T)

  # find the top preds based on function variable num_top_preds
  top_names <- names(pi_median)[order(pi_median,
                                      decreasing = T)][1:num_top_preds]
  top_names <- na.omit(top_names)

  # subset all values based on top_names
  top_pis <- ttt[ ,top_names]

  # munging and filtering for ggplot
  pi_stack <- stack(top_pis)

  # PI's have have spurious large values, NA's and NAN's
  pi_stack$values[!is.numeric(pi_stack$values)] <- NA
  pi_stack <- pi_stack[pi_stack$values < quantile(pi_stack$values,
                                                  probs = c(0.98),
                                                  na.rm = TRUE), ]
  pi_stack <- pi_stack[complete.cases(pi_stack),]

  pi_bars <- ggplot2::ggplot(pi_stack, aes(reorder(ind,
                                                   values,
                                                   FUN=median),
                                           values)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = "Relative PI", x = "")
  pi_bars
}

#' Plot PDs
#'
#' @export
plot_pds <- function(pd_name,
                     pds,
                     st_extent,
                     pointwise_pi = FALSE,
                     stixel_pds = TRUE,
                     k.cont.res = 25,
                     nnn.bs = 100,
                     equivalent.ensemble.ss = 10,
                     ci.alpha = 0.05,
                     mean.all.data = FALSE) {
  PD_MAX_RESOLUTION <- 50

  print('inside plot_pds')
  print(k.cont.res)

  # ----------------------
  # PD_NAME <- "EFFORT_HRS"
  # PD_NAME <- "UMD_FS_C4_1500_PD"
  # pipd.data.list,
  # st.extent.list = st.extent.list.E
  # pointwise.ci = T
  # k.cont.res <- 25
  # -------------------
  t.ul <- NULL
  t.ll <- NULL
  t.median <- NULL
  # Defines extent of PD predictions along independent axis
  # X.tail.level of data are trimmed the both ends where data are sparser.
  # This avoids edge effects.
  x.tail.level  <- 0.0
  # Confidence level for replicate/stixel level prediction Intervals
  # Because of the relatively small sample sizes I am using 80% PIs
  gbm.tail.prob <- 0.10  # Tail probabilty [0, 0.49]
  # Number of CV folds for data driven selection of gbm's n.trees
  # parameter. Since this doesn't work on my machine I have
  # take the simpler approach of fixing the n.trees == 200
  # This seems to work well.
  gbm.cv.folds <- 0  # >= 5 computational vs statistical efficiecy
  gbm.n.trees <- 500
  best.iter <- gbm.n.trees
  # Number of bootstrap replicasted for Conditional Mean
  #nnn.bs <- 25
  # Number of evaluation points on x-axis / indepent var
  nd.pred.size <- PD_MAX_RESOLUTION
  # ----------------------

  pd.index <- pds$centroid.date > st_extent$t.min &
              pds$centroid.date <= st_extent$t.max &
              pds$centroid.lat > st_extent$y.min &
              pds$centroid.lat <= st_extent$y.max &
              pds$centroid.lon > st_extent$x.min &
              pds$centroid.lon <= st_extent$x.max

  pd_vec <- pds[ pd.index, 	]
  # Select PD Variable
  var_pd <- pd_vec[pd_vec$V4 == pd_name,]
  # Clean
  var_pd <- var_pd[!is.na(var_pd$V5), ]
  # Each Column is one replicate estimate of PD
  # 	x = x coordinate values
  # 	y = y coordinate values
  pd.x <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.y <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.mean <- rep(NA, nrow(var_pd))
  for (rid in 1:nrow(var_pd)) {
    #rid <- 100
    pd.x[,rid] <- as.numeric(
      var_pd[rid, (PD_MAX_RESOLUTION+4):(2*PD_MAX_RESOLUTION+3)] )
    ttt <- as.numeric(var_pd[rid, 3:(PD_MAX_RESOLUTION+2)])
    pd.mean[rid] <- mean(ttt, na.rm=T)
    pd.y[,rid] <- ttt - pd.mean[rid]
  }
  pd.x <- as.data.frame(pd.x)
  pd.y <- as.data.frame(pd.y)
  # Compute Prediction Design for 1D PD
  ttt <- data.frame(
    x = stack(pd.x)[,1],
    y = stack(pd.y)[,1] )
  nd <- data.frame( x = seq(
    from = quantile(
      ttt$x, probs = x.tail.level, na.rm=T),
    to = quantile(
      ttt$x, probs = 1 - x.tail.level, na.rm=T),
    length = nd.pred.size ) )
  # PLOT STIXEL PD Replicates or just set up plot
  if (stixel_pds){
    matplot(
      jitter(as.matrix(pd.x), amount=0.00),
      pd.y,
      xlab = pd_name,
      ylab = "Deviation E(Logit Occurrence)",
      type= "l",
      #pch = 15,
      #cex = 2.0,
      lwd = 5,
      lty = 1,
      col=alpha("black", .025))
    #ylim = quantile(pd.y, probs=c(0.01, 0.99), na.rm=T))
  }
  if (!stixel_pds){
    plot(
      pd.x[,1],
      pd.y[,1],
      xlab = pd_name,
      ylab = "Deviation E(Logit Occurrence)",
      type = "n")
    #ylim = quantile(pd.y, probs=c(0.01, 0.99), na.rm=T))
  }
  abline(0,0, col="black", lwd=2)
  # -----------------
  # GBM Qunatiles
  # -----------------
  if(pointwise_pi) {
    ttt <- data.frame(
      x = stack(pd.x)[,1],
      y = stack(pd.y)[,1] )
    d.ul <- gbm::gbm(
      y ~ x,
      data = ttt,
      distribution =
        list(name="quantile",alpha=(1-gbm.tail.prob)),
      n.trees = gbm.n.trees,
      interaction.depth = 4,
      shrinkage = 0.05,
      bag.fraction = 0.5,
      train.fraction = 1.0,
      cv.folds = gbm.cv.folds,
      verbose=F,
      n.cores = 1)
    #best.iter <- gbm.perf(d.ul, method="cv", plot.it=F)
    #print(best.iter)
    d.ll <- gbm::gbm(
      y ~ x,
      data = ttt,
      distribution =
        list(name="quantile",alpha=gbm.tail.prob),
      n.trees = gbm.n.trees,
      interaction.depth = 4,
      shrinkage = 0.05,
      bag.fraction = 0.5,
      train.fraction = 1.0,
      cv.folds = gbm.cv.folds,
      verbose=F,
      n.cores = 1)
    #best.iter <- gbm.perf(d.ll, method="cv", plot.it=F)
    #cat("LL:", jjj, best.iter,"\n")
    # ------------
    #best.iter <- 200
    t.ul <- predict(d.ul,
                    newdata= nd,
                    n.trees=best.iter)
    t.ll <- predict(d.ll,
                    newdata= nd,
                    n.trees=best.iter)
    # ------------
    poly.x <- c(nd[,1], rev(nd[,1]))
    poly.y <- c(t.ll, rev(t.ul))
    polygon( poly.x, poly.y, col=alpha("red", .25), border=F)
  } # END if (pointwise.ci){

  # -----------------
  # GAM Pointwise CI for conditional mean estimate
  # via bootstrapping
  # -----------------
  if (pointwise_pi){
    # nnn.bs <- 25
    # nd.pred.size <- 50
    bs.gam.pred <- matrix(NA, nd.pred.size, nnn.bs)
    for (iii.bs in 1:nnn.bs) {
      # Take Random Sample of evaluation points
      # to account for randomness in X
      # and reduces computational time.
      # E.g. We could do this taking random
      # sample of replicates and evalution points,
      # (randomly drawing from rows and columns)
      # and set the GAM sample size at 500,
      # or the equivalent sample size when averaging
      # 10 replicate PD's each evaluated at 50 x-values.
      #
      # This suggests a nice interpretation;
      # these CI's represent the sampling variation
      # of an ensemble estimate based on a given number
      # of STIXEL PD replicates.
      #
      # equivalent.ensemble.ss <- 10
      #
      # Note that the SAME set of rows is used for every
      # column! This should be changed!
      # ------
      # Sample of STIXEL replicates
      # bs.index <- sample.int(
      # 	n = nrow(var_pd),
      # 	size = round(nrow(var_pd)*0.5),
      # 	replace = T)
      # replicate.ss <- ceiling(
      # 	equivalent.ensemble.ss*PD_MAX_RESOLUTION/length(bs.index))
      # row.index <- sample.int(
      # 	n = PD_MAX_RESOLUTION,
      # 	size = replicate.ss,
      # 	# This sets a constant fraction of available
      # 	#round(PD_MAX_RESOLUTION*0.25),
      # 	replace = F)
      random.index <-  matrix(
        (rbinom(
          n = nrow(pd.x)*ncol(pd.x),
          size = 1,
          prob = equivalent.ensemble.ss*PD_MAX_RESOLUTION/
            nrow(pd.x)/ncol(pd.x)) == 1),
        nrow(pd.x),
        ncol(pd.x))
      # dim(random.index)
      # head(random.index)
      # sum(random.index)
      ttt <- data.frame(
        x = pd.x[random.index],
        y = pd.y[random.index] )
      #x = stack(pd.x[row.index, bs.index])[,1],
      #y = stack(pd.y[row.index, bs.index])[,1] )

      s = mgcv::s
      d.gam <- mgcv::gam(
        y ~ s(x, k = k.cont.res, bs="ds", m=1),
        data = ttt,
        gamma = 1.5 )
      # summary(d.gam)
      bs.gam.pred[, iii.bs] <- predict(d.gam, newdata = nd, se=F)
      #lines(nd[,1],bs.gam.pred[, iii.bs], col=alpha("green", 0.5), lwd=2 )
    }
    t.ul <- apply(bs.gam.pred, 1, quantile, probs= 1-ci.alpha, na.rm=T)
    t.ll <- apply(bs.gam.pred, 1, quantile, probs= ci.alpha, na.rm=T)
    t.median <- apply(bs.gam.pred, 1, quantile, probs= 0.5, na.rm=T)
    # ------------
    poly.x <- c(nd[,1], rev(nd[,1]))
    poly.y <- c(t.ll, rev(t.ul))
    polygon( poly.x, poly.y, col=alpha("blue", .25), border=F)
    lines(nd[,1], t.median,
          col=alpha("darkorange", 1.0), lwd=2 )
  }
  # GAM CONDITIONAL MEAN - ALL DATA
  if (mean.all.data){
    ttt <- data.frame(
      x = stack(pd.x)[,1],
      y = stack(pd.y)[,1] )
    d.gam <- mgcv::gam(
      y ~ s(x, k = k.cont.res, bs="ds", m=1),
      data = ttt,
      gamma = 1.5 )
    # summary(d.gam)
    p.gam <- predict(d.gam,
                     newdata = nd,
                     se=T)
    polygon(
      x = c(nd[,1], rev(nd[,1])),
      y = c( p.gam$fit + 2*p.gam$se.fit,
             rev(p.gam$fit - 2*p.gam$se.fit) ),
      col = alpha("lightblue", 0.25),
      border = NA)
    lines( nd[,1], p.gam$fit, col = alpha("yellow", 0.75), lwd=3)
  }
  return( list(
    t.median = t.median,
    t.ul = t.ul,
    t.ll = t.ll))
}
