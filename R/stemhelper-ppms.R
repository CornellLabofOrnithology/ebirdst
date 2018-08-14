#' Load Predictive Performance Metrics Data
#'
#' Internal function that loads test data and filters it against ensemble
#' support data.
#'
#' @param path character; Full path to single species STEM results.
#'
#' @return list with data.frame of ppm_data and data.frame of interpolated
#' daily ensemble support values for the Western Hemisphere.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' ppm_data <- load_ppm_data(sp_path)
#' }
load_ppm_data <- function(path) {
  # load the test data and assign names
  test_file <- paste(path,
                     "/results/abund_preds/unpeeled_folds/test.pred.ave.txt",
                     sep = "")

  if(!file.exists(test_file)) {
    stop("*_erd.test.data.csv file does not exist in the /data directory.")
  }

  ppm_data <- data.table::fread(test_file, showProgress = FALSE)
  ppm_names <- c("data.type", "row.id", "lon", "lat", "date", "obs", "pi.mean",
                 "pi.90", "pi.10", "pi.se", "pi.mu.mean", "pi.mu.90",
                 "pi.mu.10", "pi.mu.se", "pat", "pi.es")
  names(ppm_data) <- ppm_names

  # load ensemble support values that define weekly extent of analysis
  es_dir <- paste(path,
                  "/results/tifs/presentation/",
                  "abundance_ensemble_support_values/",
                  sep = "")
  es_files <- list.files(es_dir)

  eoa_es_data <- data.frame(Week = NA,
                            WesternHemisphere = NA,
                            Pat = NA)

  # TODO...how does this behave if there are less than 52 weeks of data?
  for (iii in 1:length(es_files)) {
    ttt <- utils::read.csv(paste(es_dir, es_files[iii], sep = ""))
    eoa_es_data[iii, 1] <- as.character(ttt[1, 2])
    eoa_es_data[iii, 2] <- ttt[1, 3]
    eoa_es_data[iii, 3] <- ttt[1, 4]
  }

  # add day of year
  eoa_es_data$DOY <- as.numeric(format(strptime(x = eoa_es_data$Week,
                                                format = "%m-%d"),
                                       "%j"))

  # smooth the extent of estimate ensemble support values to day

  # init
  pred_DOY <- data.frame(DOY = c(1:366))
  eoa_WesternHemisphere <- rep(50, nrow(pred_DOY))
  eoa_Pat <- rep(NA, nrow(pred_DOY))

  slope = (25 - 3)/(52 - 10)

  # check for NA data
  na_total <- sum(!is.na(eoa_es_data$Pat))

  if(na_total > 1) {
    # Treat missing Pat values as the minimum support value
    pat_na_replace <- min(eoa_es_data$Pat, na.rm = TRUE)
    eoa_es_data$Pat[is.na(eoa_es_data$Pat)] <- pat_na_replace

    if(na_total > 9) {
      # GAM Smooth EOA ES values down to Daily
      y = -((52 * slope - na_total * slope) - 25)

      s = mgcv::s
      na_model <- mgcv::gam(Pat ~ s(DOY,
                                             k = round(y),
                                             bs = "cp",
                                             m = 1),
                            gamma = 1.5,
                            data = eoa_es_data,
                            knots = list(DOY = c(1,366)))
    } else {
      # do a loess
      na_model <- stats::loess(Pat ~ DOY, eoa_es_data)
    }

    eoa_Pat <- stats::predict(na_model, newdata = pred_DOY)
  }

  # package and return
  eoa_es_daily <- data.frame(DOY = pred_DOY,
                             WesternHemisphere = eoa_WesternHemisphere,
                             Pat = eoa_Pat)

  return(list(ppm_data = ppm_data,
              eoa_es_daily = eoa_es_daily))
}

#' Computes the Predictive Performance Metrics for a spatiotemporal extent
#'
#' Loads test data and ensemble support values and then calculates the
#' predictive performance metrics within a spatiotemporal extent defined by
#' `st_extent`. Use this function directly to access the computed metrics, or
#' use `plot_all_ppms()` or `plot_binary_by_time()` to summarize the metrics.
#'
#' @param path character; Full path to single species STEM results.
#' @param st_extent list; st_extent list for spatiotemporal filtering.
#'
#' @return list of data.frames: binary_stats, occ_stats, count_stats
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' ppms <- compute_ppms(path = sp_path, st_extent = ne_extent)
#' }
compute_ppms <- function(path, st_extent = NA) {

  poisson.dev <- function(obs, pred) {
    dev.mean <- NA
    dev.model <- NA
    dev.explained <- NA

    meanpred <- mean(obs, na.rm = TRUE)

    dev.mean <- 2 * sum(obs * log(ifelse( obs == 0, 1, obs / meanpred)) -
                          (obs - meanpred),
                        na.rm = TRUE)
    dev.model <- 2 * sum(obs * log(ifelse( obs == 0, 1, obs / pred)) -
                           (obs - pred),
                         na.rm = TRUE)
    dev.explained <- 1 - dev.model / dev.mean

    return(c(dev.model, dev.mean, dev.explained))
  }

  bernoulli.dev <- function(obs, pred) {
    dev.mean <- NA
    dev.model <- NA
    dev.explained <- NA

    meanpred <- mean(obs, na.rm = TRUE)

    dev.mean <- 2 * sum(obs * log(ifelse(meanpred == 0, 1, meanpred)) +
                          (1 - obs) *
                          log(ifelse(1 - meanpred == 0, 1, 1 - meanpred)),
                        na.rm = TRUE)
    dev.model <- 2 * sum(obs * log(ifelse(pred == 0, 1, pred)) +
                           (1 - obs) *
                           log(ifelse(1 - pred == 0, 1, 1 - pred)),
                         na.rm = TRUE)
    dev.explained <- 1 - dev.model / dev.mean

    return(c(dev.model, dev.mean, dev.explained))
  }

  balanced.binary.sample <- function(
    binary.outcome,
    x,
    y,
    xlim = NA,
    ylim = NA,
    nx = 50,
    ny = 50,
    min.class = 0.25,
    neg.size = NA) {

    sample.grid.cell <- function(
      xxx,
      yyy,
      xlim = c(NA,NA),
      ylim = c(NA,NA),
      nx = 64,
      ny = 64,
      jitter = F,
      size = 1,
      replace = F ){

      lookup.grid.cell <- function(
        xxx,
        yyy,
        xlim = c(NA,NA),
        ylim = c(NA,NA),
        nx = 64,
        ny = 64,
        jitter = F ){
        # -----------------------------------
        cell.number <- rep(NA, length(xxx))
        x.ll <- y.ll <- NA
        xxx.width <- yyy.width <- NA
        if (length(xxx)==length(yyy) & length(unique(xxx))>1){
          if (any(is.na(xlim))) xlim <- range(xxx, na.rm=T)
          if (any(is.na(ylim))) ylim <- range(yyy, na.rm=T)
          xxx.width <- abs(xlim[2]-xlim[1])/nx
          yyy.width <- abs(ylim[2]-ylim[1])/ny
          # Lower Left Corner
          x.ll <- min(xlim)
          y.ll <- min(ylim)
          # If Jittered, domain is bigger & number of grid cells increases
          if (jitter){
            x.ll <- x.ll - stats::runif(1)*xxx.width
            y.ll <- y.ll - stats::runif(1)*yyy.width
            nx <- 1 + (max(xlim) - min( c(x.ll, xlim) )) %/% xxx.width
            ny <- 1 + (max(ylim) - min( c(y.ll, ylim) )) %/% yyy.width
          }
          # ID data within grid
          ingrid.index <-
            xxx >= x.ll & xxx <= x.ll + nx*xxx.width &
            yyy >= y.ll & yyy <= y.ll + ny*yyy.width
          if (sum(ingrid.index)>0){
            # Assign Row, Column, Grid Cell Number
            col.number <- 1 + ( xxx[ingrid.index] - x.ll ) %/% xxx.width
            row.number <- 1 + ( yyy[ingrid.index] - y.ll ) %/% yyy.width
            cell.number[ingrid.index] <- col.number + (row.number-1)*nx
          }
        }
        return( list(
          cell.number=cell.number,
          # jittered bounding box
          bb = matrix(
            c(x.ll, y.ll, x.ll + nx*xxx.width, y.ll + ny*yyy.width), 2, 2,
            byrow=T, dimnames=list(c("ll", "ur"), c("xxx", "yyy"))),
          nx = nx,
          ny = ny,
          xwidth = xxx.width,
          ywidth = yyy.width  ))
      }

      # Stratified sample over Grid Cell Number
      sample_fun <- function(x, size, replace){
        # Cells without samples are excluded in the
        # tapply call - if (length(x)==0) return(NA)
        # Cells with a single sample cause problems, see help(sample)
        # So, I am going to handle this situation "by hand"
        result <- rep(NA, size)
        if (length(x)==1 & replace==F) {
          #cat("sf: length(x)==1 & replace==F",x,"\n")
          result <- rep(NA, size)
          result[1] <- x
        }
        if (length(x)==1 & replace==T) {
          #cat("sf: length(x)==1 & replace==T",x,"\n")
          result <- rep(x, size)
        }
        if (length(x)>1 & replace == F & size > length(x) ){
          result <- rep(NA, size)
          result[1:length(x)] <- x
        }
        if (length(x)>1 & replace == F & size <= length(x) ){
          result <- sample(x=x, size=size, replace=replace)
        }
        if (length(x)>1 & replace == T ){
          result <- sample(x=x, size=size, replace=replace)
        }
        return(result)
      }
      lgc <- lookup.grid.cell(
        xxx, yyy, xlim, ylim, nx, ny, jitter)
      n.index <- tapply(
        c(1:length(xxx))[!is.na(lgc$cell.number)],
        as.factor(lgc$cell.number[!is.na(lgc$cell.number)]),
        sample_fun, size, replace)
      n.index <- plyr::rbind.fill.matrix(n.index)
      return(list(
        cell.number = lgc$cell.number,
        sample.index = n.index,
        bb = lgc$bb,
        nx = lgc$nx,
        ny = lgc$ny,
        xwidth = lgc$xwidth,
        ywidth = lgc$ywidth  ))
    }

    # ----------------------------------------------------------------------
    sgc.pos.nindex <- NULL
    sgc.neg.nindex <- NULL
    jjj.max <- 20 # Maximum spatially balanced oversampling iterations
    if (is.na(neg.size)) neg.size <- 1
    if (any(is.na(xlim))) xlim = range(x, na.rm=T)
    if (any(is.na(ylim))) ylim = range(y, na.rm=T)
    if (length(binary.outcome)==length(x) &
        length(x)==length(y) &
        length(binary.outcome)>1){
      # -----------------------
      # Check for negatives
      if (sum(binary.outcome==0) > 0) {
        ttt.nindex <- c(1:length(binary.outcome))[ binary.outcome == 0 ]
        xxx <- x[ttt.nindex]
        yyy <- y[ttt.nindex]
        sgc.neg <- sample.grid.cell(
          xxx, yyy, xlim, ylim, nx, ny, jitter = T,
          size = neg.size,
          replace = F )
        # Remember to strip out NA's from grid
        sgc.neg$sample.index <- sgc.neg$sample.index[!is.na(sgc.neg$sample.index)]
        # Index back to original data.frame
        sgc.neg.nindex <- ttt.nindex[ sgc.neg$sample.index ]
      }
      # -----------------------
      # Check for positives
      if (sum(binary.outcome>0) > 0) {
        ttt.nindex <- c(1:length(binary.outcome))[ binary.outcome != 0 ]
        xxx <- x[ttt.nindex]
        yyy <- y[ttt.nindex]
        # Check Positive Classs Proportion
        # Compared to Balanced Neg Sample
        ncut <- round(min.class/(1-min.class)*length(sgc.neg.nindex))
        if ( length(ttt.nindex) < ncut) {
          # Oversample positives
          sgc.pos.nindex <-  sample(
            ttt.nindex,
            size = 1+ncut,
            replace = T)
        }
        # Else Spatial Sample
        if ( length(ttt.nindex) >= ncut ) {
          sgc.pos <- sample.grid.cell(
            xxx, yyy, xlim, ylim, nx, ny, jitter = T, size = 1, replace = F )
          sgc.pos$sample.index <- sgc.pos$sample.index[!is.na(sgc.pos$sample.index)]
          sgc.pos.nindex <- ttt.nindex[ sgc.pos$sample.index ]
          # Check proportion
          # If less than min.class then increase sample per cell until
          # minimum proportion achieved.
          jjj <- 1
          while ( length(sgc.pos.nindex) < ncut & jjj < jjj.max ) {
            jjj <- jjj+1
            #print(jjj)
            sgc.pos <- sample.grid.cell(
              xxx, yyy, xlim, ylim, nx, ny, jitter = T,
              size = jjj,
              replace = F )
            sgc.pos$sample.index <- sgc.pos$sample.index[!is.na(sgc.pos$sample.index)]
            sgc.pos.nindex <- ttt.nindex[ sgc.pos$sample.index ]
          }
        } # END else spatial sample
      } # End check for positives
    }# Check input parameter lengths
    return(
      list(
        pos.nindex = sgc.pos.nindex,
        neg.nindex = sgc.neg.nindex ))
  }

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  e <- load_config(path)

  # load the test data and assign names
  test_file <- paste(path,
                     "/results/abund_preds/unpeeled_folds/test.pred.ave.txt",
                     sep = "")

  if(!file.exists(test_file)) {
    stop("*_erd.test.data.csv file does not exist in the /data directory.")
  }

  ppm_data <- data.table::fread(test_file, showProgress = FALSE)
  ppm_names <- c("data.type", "row.id", "lon", "lat", "date", "obs", "pi.mean",
                 "pi.90", "pi.10", "pi.se", "pi.mu.mean", "pi.mu.90",
                 "pi.mu.10", "pi.mu.se", "pat", "pi.es")
  names(ppm_data) <- ppm_names

  # static vars
  n_mc <- 25
  occ_binary <- TRUE
  occ_prob <- TRUE
  count <- TRUE
  # Min Occ SS within range
  occ_min_ss <- 50
  occ_min_mean <- 0.01
  # Min Occ SS within range
  count_min_ss <- 50
  count_min_mean <- 0.25
  # Set Count Max
  count_max <- NA

  # names -------------------
  # Binary / Range PPMs
  binary_stat_names <- c("mc", "ss", "mean", "AUC", "PCC", "Kappa", "B.DE",
                         "Sensitivity", "Specificity")
  # Within range: Occurrence Rate PPMs
  occ_stat_names <- c("mc", "ss", "mean", "threshold", "AUC", "PCC", "Kappa",
                      "B.DE", "Sensitivity", "Specificity")
  # Within range: Expected Count PPMs
  count_stat_names <- c("mc", "ss", "mean", "P.DE.abund", "P.DE.occ",
                        "Spearman.abund", "Spearman.occ")

  # ------------------------------------------------------------------
  # Compute MC Sample of PPMs for ST Subsets
  # ------------------------------------------------------------------
  # Extract ST Subset

  if(!all(is.na(st_extent))) {
    st_data <- st_extent_subset(ppm_data, st_extent)
  } else {
    st_data <- ppm_data
  }

  rm(ppm_data)

  # stack stem for the week
  time_stack <- stack_stem(path, variable = "abundance_umean", res = "high",
                           year = 2016, st_extent = st_extent,
                           add_zeroes = TRUE)

  # TODO add logic for averaging if there are more than one

  # spatialize the st_data
  st_data_sp <- sp::SpatialPointsDataFrame(st_data[, c("lon", "lat")],
                                           st_data,
                                           proj4string =
                                             sp::CRS("+init=epsg:4326"))
  sinu <- sp::CRS(sp::proj4string(stemhelper::template_raster))
  st_data_prj <- sp::spTransform(st_data_sp, sinu)
  rm(st_data_sp)

  st_data_e <- raster::extract(time_stack, st_data_prj)

  print(sum(is.na(st_data_e)))

  st_data <- st_data[!is.na(st_data_e), ]
  rm(st_data_e, time_stack, st_data_prj, sinu)

  # Split data into within extent and out of extent
  st_data_zeroes <- st_data[st_data$pi.es < 75, ]
  st_data <- st_data[st_data$pi.es >= 75, ]

  # FDR
  binom_test_p <- function(x) {
    binom.test(round(as.numeric(x["pat"]) * as.numeric(x["pi.es"]), 0),
               as.numeric(x["pi.es"]),
               0.10,
               alternative = "greater")$p.value
  }

  p_values <- apply(st_data, 1, binom_test_p)
  p_adj <- p.adjust(p_values, "fdr")

  # Add Binary Prediction
  st_data$binary <- as.numeric(p_adj < 0.001)

  # Readd test data out of extent with binary as 0
  st_data_zeroes$binary <- 0

  st_data <- rbind(st_data, st_data_zeroes)

  # Remove rows where binary = NA (inherited from PAT = NA)
  st_data <- st_data[!is.na(st_data$binary), ]

  if(nrow(st_data) == 0) {
    return(NA)
  }

  # -------------------------------------------------------------
  # Monte Carlo Samples for PPMs Binary, Spatially Balanced Sample
  # (NO!!! oversampling +'s')
  # -------------------------------------------------------------
  binary_stats <- matrix(NA, n_mc, length(binary_stat_names))
  binary_stats <- as.data.frame(binary_stats)
  names(binary_stats) <- binary_stat_names

  occ_stats <- matrix(NA, n_mc, length(occ_stat_names))
  occ_stats <- as.data.frame(occ_stats)
  names(occ_stats) <- occ_stat_names

  count_stats <- matrix(NA, n_mc, length(count_stat_names))
  count_stats <- as.data.frame(count_stats)
  names(count_stats) <- count_stat_names

  for(iii.mc in 1:n_mc) {
    # -------------------------------------------------------------
    # Sampling
    # -------------------------------------------------------------
    # Spatially Balanced Case Control Sampling
    # (NOTE NO!!! oversampling)
    #bbs <- balanced.binary.sample(binary.outcome = as.numeric(st_data$obs > 0),
    #                              x = st_data$lon,
    #                              y = st_data$lat,
    #                              xlim = range(st_data$lon, na.rm = TRUE),
    #                              ylim = range(st_data$lat, na.rm = TRUE),
    #                              nx = 50,
    #                              ny = 50,
    #                              min.class = 0.0,
    #                              neg.size = 1)

    # projected bbs
    st_data_sp <- sp::SpatialPointsDataFrame(st_data[, c("lon", "lat")],
                                             st_data,
                                             proj4string =
                                               sp::CRS("+init=epsg:4326"))
    sinu <- sp::CRS(sp::proj4string(stemhelper::template_raster))
    st_data_prj <- sp::spTransform(st_data_sp, sinu)

    xrange <- range(st_data_prj@coords[, 1], na.rm = TRUE)
    yrange <- range(st_data_prj@coords[, 2], na.rm = TRUE)

    xwidth <- xrange[2] - xrange[1]
    yheight <- yrange[2] - yrange[1]

    bbs <- balanced.binary.sample(binary.outcome = as.numeric(st_data$obs > 0),
                                  x = st_data_prj@coords[, 1],
                                  y = st_data_prj@coords[, 2],
                                  xlim = xrange,
                                  ylim = yrange,
                                  nx = xwidth/10000,
                                  ny = yheight/10000,
                                  min.class = 0.0,
                                  neg.size = 1)

    # Index back to full vector
    sample.nindex <- c(1:nrow(st_data))[c(bbs$pos.nindex, bbs$neg.nindex)]

    print(length(bbs$pos.nindex))
    print(length(bbs$neg.nindex))
    print((length(bbs$pos.nindex) / length(bbs$neg.nindex)) * 100)

    ttt.data <- st_data[sample.nindex, ]

    # -------------------------------------------------------------
    # Binary Occupancy PPMs
    # -------------------------------------------------------------
    if(occ_binary) {
      binary_stats$ss[iii.mc] <- nrow(ttt.data)
      binary_stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs > 0))

      if(nrow(ttt.data) >= occ_min_ss &
         mean(as.numeric(ttt.data$obs > 0)) >= occ_min_mean) {
        # ------------
        pa.df <- data.frame(
          "space.holder.tag",
          obs = as.numeric(ttt.data$obs > 0),
          ppp = ttt.data$binary)
        pa.cmx <- PresenceAbsence::cmx(pa.df, na.rm = T)

        binary_stats$mc[iii.mc] <- iii.mc
        pa_sensitivity <- PresenceAbsence::sensitivity(pa.cmx, st.dev = FALSE)
        binary_stats$Sensitivity[iii.mc] <- pa_sensitivity
        pa_specificity <- PresenceAbsence::specificity(pa.cmx, st.dev = FALSE)
        binary_stats$Specificity[iii.mc] <- pa_specificity
        pa_kappa <- PresenceAbsence::Kappa(pa.cmx, st.dev = FALSE)
        binary_stats$Kappa[iii.mc] <- pa_kappa
        pa_pcc <- PresenceAbsence::pcc(pa.cmx, st.dev = FALSE)
        binary_stats$PCC[iii.mc] <- pa_pcc
        pa_auc <- PresenceAbsence::auc(pa.df, na.rm = TRUE, st.dev = FALSE)
        binary_stats$AUC[iii.mc] <- pa_auc
        bde <- as.numeric(bernoulli.dev(obs = pa.df$obs, pred = pa.df$ppp))[3]
        binary_stats$B.DE[iii.mc] <- bde
      } # END MIN SS
    } # occ.binary

    # -------------------------------------------------------------
    # Within range, Occupancy Rate PPMs
    # -------------------------------------------------------------
    if(occ_prob | count) {
      # Limit ttt.data to within Range
      ttt.data <- ttt.data[ttt.data$binary > 0, ]
    }

    if(occ_prob) {
      # Remove Rows with missing pi.mean
      ttt.data <- ttt.data[!is.na(ttt.data$pi.mean), ]
      ttt.data <- ttt.data[!is.nan(ttt.data$pi.mean), ]

      # Min Occ SS
      occ_stats$ss[iii.mc] <- nrow(ttt.data)
      occ_stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs > 0))

      if(nrow(ttt.data) >= occ_min_ss &
          mean(as.numeric(ttt.data$obs > 0)) >= occ_min_mean ) {
        # ------------
        pa2.df <- data.frame("space.holder.tag",
                             obs = as.numeric(ttt.data$obs > 0),
                             ppp = ttt.data$pi.mean)
        pa.metrics <- PresenceAbsence::presence.absence.accuracy(pa2.df,
                                                                 threshold = 0.5,
                                                                 na.rm = TRUE,
                                                                 st.dev = FALSE)
        # -----------------
        # Note, we use fixed estimate of the threshold instead of the
        # data driven method,
        #	optimal_thresh_position <- which.max(pa.metrics$Kappa)
        # We do this to preserve the indepednece of test set.
        # It makes a little difference, but not much.
        # -----------------
        optimal_thresh_position <- 1
        # ------------------
        occ_stats$mc[iii.mc] <- iii.mc
        occ_thresh <- as.numeric(pa.metrics$threshold[optimal_thresh_position])
        occ_stats$threshold[iii.mc] <- occ_thresh
        occ_pcc <- as.numeric(pa.metrics$PCC[optimal_thresh_position])
        occ_stats$PCC[iii.mc] <- occ_pcc
        occ_sens <- as.numeric(pa.metrics$sensitivity[optimal_thresh_position])
        occ_stats$Sensitivity[iii.mc] <- occ_sens
        occ_spec <- as.numeric(pa.metrics$specificity[optimal_thresh_position])
        occ_stats$Specificity[iii.mc] <- occ_spec
        occ_kappa <- as.numeric(pa.metrics$Kappa[optimal_thresh_position])
        occ_stats$Kappa[iii.mc] <- occ_kappa
        occ_auc <- as.numeric(pa.metrics$AUC[optimal_thresh_position])
        occ_stats$AUC[iii.mc] <- occ_auc
        occ_bde <- as.numeric(bernoulli.dev(obs = pa2.df$obs,
                                            pred = pa2.df$ppp))[3]
        occ_stats$B.DE[iii.mc] <- occ_bde
      } # END Min Occ SS
    } # END occ.prob
    # -------------------------------------------------------------
    # Within range, Expected Count PPMs
    # -------------------------------------------------------------
    if(count) {
      # Remove Rows with missing pi.mu.mean
      ttt.data <- ttt.data[!is.na(ttt.data$pi.mu.mean), ]
      ttt.data <- ttt.data[!is.nan(ttt.data$pi.mu.mean), ]

      # Count Max
      if(!is.na(count_max)) {
        ttt.data$obs[ttt.data$obs > count_max] <- count_max
      }

      # Min Count SS
      count_stats$ss[iii.mc] <- nrow(ttt.data)
      count_stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs))

      # --------------
      if(nrow(ttt.data) >= count_min_ss &
         mean(as.numeric(ttt.data$obs)) >= count_min_mean ) {
        # --------------
        count_stats$mc[iii.mc] <- iii.mc
        # Poisson Deviance
        pdev <- as.numeric(poisson.dev(obs = ttt.data$obs,
                                       pred = ttt.data$pi.mu.mean))
        count_stats$P.DE.abund[iii.mc] <- pdev[3]

        pdev <- as.numeric(poisson.dev(obs = ttt.data$obs,
                                       pred = ttt.data$pi.mean))
        count_stats$P.DE.occ[iii.mc] <- pdev[3]

        # Spearman's Rank Correlations
        count_stats$Spearman.abund[iii.mc] <- stats::cor(ttt.data$pi.mu.mean,
                                                         ttt.data$obs,
                                                         method = "spearman")
        count_stats$Spearman.occ[iii.mc] <- stats::cor(ttt.data$pi.mean,
                                                       ttt.data$obs,
                                                       method = "spearman")
      } # END Min Count SS
    } # END count
  } # END iii.mc

  return(list(binary_stats = binary_stats,
              occ_stats = occ_stats,
              count_stats = count_stats))
}

#' Plot binary occurrence metrics by time
#'
#' For a specified number of time periods (ideally weeks or months), plots one
#' of four (Kappa, AUC, Sensitivity, Specificity) box plots. Use default,
#' without `st_extent` to see range-wide predictive performance, or provide
#' an `st_extent` to see performance within spatiotemporal extent.
#'
#' @param path character; Full path to single species STEM results.
#' @param metric character; One of: "Kappa", "AUC", "Sensitivity",
#' "Specificity".
#' @param n_time_periods int; Greater than 1 (e.g., 52 for weeks).
#' @param st_extent list; st_extent list for spatiotemporal filtering.
#'
#' @return Plot of metric box plots by time.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#'
#' plot_binary_by_time(path = sp_path, metric = "Kappa", n_time_periods = 12)
#' }
plot_binary_by_time <- function(path,
                                metric,
                                n_time_periods = 52,
                                st_extent = NA) {

  if(!(metric %in% c("Kappa", "AUC", "Sensitivity", "Specificity"))) {
    stop(paste("Predictive performance metric must be one of: ",
               "Kappa, AUC, Sensitivity, Specificity.",
               sep = ""))
  }

  if(n_time_periods < 2) {
    stop("n_time_periods argument must be more than 1")
  }

  seasonal_ppms <- list(NA)

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  if(all(is.na(st_extent))) {
    st_extent <- list()
    st_extent$type = "rectangle"
    st_extent$lat.min <- -90
    st_extent$lat.max <- 90
    st_extent$lon.min <- -180
    st_extent$lon.max <- 180
  }

  for(i_t in 1:n_time_periods) {
    st_extent$t.min <- (i_t-1)/n_time_periods
    st_extent$t.max <- i_t/n_time_periods

    s_ppms <- compute_ppms(path, st_extent)
    seasonal_ppms[[i_t]] <- s_ppms
  }

  ttt <- NULL

  for (iii.time in 1:n_time_periods) {
    if(!all(is.na(seasonal_ppms[[iii.time]]))) {
      values <- seasonal_ppms[[iii.time]][[1]][, c(metric)]
      if(sum(is.na(values)) == length(values)) {
        values <- 0
      }

      # Add mean resp
      ttt <- rbind(ttt,
                   cbind(rep(iii.time,
                             length(values)),
                         values,
                         seasonal_ppms[[iii.time]][[1]][, 3]))
    }
  }

  ttt <- as.data.frame(ttt)
  names(ttt) <- c("Week", "Values", "mean.resp")

  ttt$Date <- apply(ttt, 1, FUN = function(x) {
    strftime(as.Date((x["Week"] / n_time_periods) * 366,
                     origin = as.Date('2013-01-01')),
             format = "%Y-%m-%d")
  })

  ttt$Date <- as.Date(ttt$Date)



  # -----------------
  # Plot
  # -----------------

  if(metric == "AUC") {
    metric_ylim <- c(0.5, 1)
  } else {
    metric_ylim <- c(0, 1)
  }

  bp <- ggplot2::ggplot(ttt,
                        ggplot2::aes(ttt$Date, ttt$Values, group = ttt$Date)) +
    ggplot2::geom_boxplot() +
    ggplot2::ylim(metric_ylim) +
    ggplot2::xlab("Date") +
    ggplot2::ggtitle(metric) +
    ggplot2::scale_x_date(date_labels = "%b",
                          limits = c(as.Date("2013-01-01"),
                                     as.Date("2013-12-31")),
                          date_breaks = "1 month") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  bp

}

#' Plot all predictive performance metrics
#'
#' For a spatiotemporal extent, plots bar plots for all available predictive
#' performance metrics within three sets: Binary Occurrence, Occurrence
#' Probability, and Abundance.
#'
#' @param path character; Full path to single species STEM results.
#' @param st_extent list; st_extent list for spatiotemporal filtering.
#'
#' @return Plot of metric box plots by category
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#'
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' plot_all_ppms(path = sp_path, st_extent = ne_extent)
#' }
plot_all_ppms <- function(path, st_extent) {
  if(all(is.na(st_extent))) {
    stop("Must provide a complete spatiotemporal extent.")
  } else {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  if(is.null(st_extent$t.min) | is.null(st_extent$t.max)) {
    stop("Must provide temporal limits (t.min, t.max) in st_extent.")
  }

  ppm_data <- compute_ppms(path, st_extent)

  ppm_data$binary_stats$type <- "Occ Binary"
  ppm_data$occ_stats$type <- "Occ Probability"
  ppm_data$count_stats$type <- "Abundance"

  all_ppms <- dplyr::bind_rows(ppm_data)

  # columns: metric, values, type
  all_ppms_melt <- reshape2::melt(all_ppms, id = c("type"))
  all_ppms_melt$variable <- as.character(all_ppms_melt$variable)
  all_ppms_melt <- all_ppms_melt[!(all_ppms_melt$variable %in%
                                     c("mc", "ss", "mean", "threshold",
                                       "P.DE.occ", "Spearman.occ")), ]
  all_ppms_melt <- all_ppms_melt[!(all_ppms_melt$type == "Occ Binary" &
                                     all_ppms_melt$variable == "B.DE"), ]
  all_ppms_melt <- all_ppms_melt[!is.na(all_ppms_melt$value), ]

  # Separate
  all_ppms_melt_ob <- all_ppms_melt[all_ppms_melt$type == "Occ Binary", ]

  obp <- ggplot2::ggplot(all_ppms_melt_ob,
                        ggplot2::aes(all_ppms_melt_ob$variable,
                                     all_ppms_melt_ob$value,
                                     group = all_ppms_melt_ob$variable)) +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab("Occ Binary") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.y = ggplot2::element_blank())
  obp

  all_ppms_melt_op <- all_ppms_melt[all_ppms_melt$type == "Occ Probability", ]
  all_ppms_melt_op$variable <- factor(all_ppms_melt_op$variable,
                                      levels = c("AUC", "Kappa", "B.DE", "PCC",
                                                 "Sensitivity", "Specificity"))
  # negative check of B.DE
  medbde <- stats::median(all_ppms_melt_op[
    all_ppms_melt_op$variable == "B.DE", ]$value)

  bderep <- data.frame(type = "Occ Probability",
                       variable = "B.DE",
                       value = 0,
                       stringsAsFactors = FALSE)

  opp <- ggplot2::ggplot(all_ppms_melt_op,
                         ggplot2::aes(all_ppms_melt_op$variable,
                                      all_ppms_melt_op$value,
                                      group = all_ppms_melt_op$variable)) +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab("Occ Probability") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.y = ggplot2::element_blank())

  # if the median B.DE value is below 0, put a red X
  if(medbde < 0) {
    opp <- opp + ggplot2::geom_point(data = bderep,
                                     ggplot2::aes(bderep$variable,
                                                  bderep$value,
                                                  group = bderep$variable),
                                     shape = 4,
                                     size = 10,
                                     color = 'red')
  }

  opp

  all_ppms_melt_ab <- all_ppms_melt[all_ppms_melt$type == "Abundance", ]
  all_ppms_melt_ab[all_ppms_melt_ab$variable == "Spearman.abund", ]$variable <- "Spearman"

  # negative check of P.DE
  medpde <- stats::median(all_ppms_melt_ab[
    all_ppms_melt_ab$variable == "P.DE.abund", ]$value)

  pderep <- data.frame(type = "Abundance",
                       variable = "P.DE.abund",
                       value = 0,
                       stringsAsFactors = FALSE)

  abp <- ggplot2::ggplot(all_ppms_melt_ab,
                         ggplot2::aes(all_ppms_melt_ab$variable,
                                      all_ppms_melt_ab$value,
                                      group = all_ppms_melt_ab$variable)) +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab("Abundance") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       size = 8),
                   axis.title.y = ggplot2::element_blank())

  # if the median B.DE value is below 0, put a red X
  if(medpde < 0) {
    abp <- abp + ggplot2::geom_point(data = pderep,
                                     ggplot2::aes(pderep$variable,
                                                  pderep$value,
                                                  group = pderep$variable),
                                     shape = 4,
                                     size = 10,
                                     color = 'red')
  }

  abp

  three_plots <- list(obp, opp, abp)
  gridExtra::grid.arrange(grobs = three_plots, ncol = 3)
}
