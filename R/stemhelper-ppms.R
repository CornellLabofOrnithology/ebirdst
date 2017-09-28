#' Internal lookup.grid.cell function
#'
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
      x.ll <- x.ll - runif(1)*xxx.width
      y.ll <- y.ll - runif(1)*yyy.width
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

#' Internal sample.grid.cell function
#'
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

  # Stratified sample over Grid Cell Number
  sample_fun <- function(x, size, replace){
    # Cells without samples are excluded in the tapply call - if (length(x)==0) return(NA)
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

#' Internal balanced binary sample function
#'
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

#' Load PPM data
#'
#' @export
load_ppm_data <- function(path) {
  # load the test data and assign names
  test_file <- paste(path,
                     "/results/abund_preds/unpeeled_folds/test.pred.ave.txt",
                     sep = "")

  if(!file.exists(test_file)) {
    stop("*_erd.test.data.csv file does not exist in the /data directory.")
  }

  ppm_data <- data.table::fread(test_file)
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
                            NorthAmerica = NA,
                            SouthAmerica = NA)

  for (iii in 1:length(es_files)) {
    ttt <- read.csv(paste(es_dir, es_files[iii], sep = ""))
    eoa_es_data[iii, 1] <- as.character(ttt[1, 2])
    eoa_es_data[iii, 2:3] <- ttt[1, 3:4]
  }

  # add day of year
  eoa_es_data$DOY <- as.numeric(format(strptime(x = eoa_es_data$Week,
                                                format = "%m-%d"),
                                       "%j"))

  # smooth the extent of estimate ensemble support values to day

  # init
  pred_DOY <- data.frame(DOY = c(1:366))
  eoa_NorthAmerica <- rep(NA, nrow(pred_DOY))
  eoa_SouthAmerica <- rep(NA, nrow(pred_DOY))

  # check for NA data
  if(sum(!is.na(eoa_es_data$NorthAmerica)) > 10) {
    # Treat missing ES values as the minimum support value
    na_na_replace <- min(eoa_es_data$NorthAmerica, na.rm = TRUE)
    eoa_es_data$NorthAmerica[is.na(eoa_es_data$NorthAmerica)] <- na_na_replace

    # GAM Smooth EOA ES values down to Daily
    s = mgcv::s
    na_gam <- mgcv::gam(NorthAmerica ~ s(DOY, k = 25, bs = "cp", m = 1),
                        gamma = 1.5,
                        data = eoa_es_data,
                        knots = list(DOY = c(1,366)))
    eoa_NorthAmerica <- predict(na_gam, newdata = pred_DOY)
  }

  if(sum(!is.na(eoa_es_data$SouthAmerica)) > 10) {
    # Treat missing ES values as the minimum support value
    sa_na_replace <- min(eoa_es_data$SouthAmerica, na.rm = TRUE)
    eoa_es_data$SouthAmerica[is.na(eoa_es_data$SouthAmerica)] <- sa_na_replace

    # GAM Smooth EOA ES values down to Daily
    s = mgcv::s
    sa_gam <- mgcv::gam(SouthAmerica ~ s(DOY, k = 25, bs = "cp", m = 1),
                        gamma = 1.5,
                        data = eoa_es_data,
                        knots = list(DOY = c(1,366)))
    eoa_SouthAmerica <- predict(sa_gam, newdata = pred_DOY)
  }

  # package and return
  eoa_es_daily <- data.frame(DOY = pred_DOY,
                             NorthAmerica = eoa_NorthAmerica,
                             SouthAmerica = eoa_SouthAmerica  )

  return(list(ppm_data = ppm_data,
              eoa_es_daily = eoa_es_daily))
}

#' Internal poisson.dev function
#'
poisson.dev <- function(obs, pred){
  dev.mean <- NA
  dev.model <- NA
  dev.explained <- NA
  meanpred <- mean(obs,na.rm=T)
  dev.mean <- 2*sum( obs *
                       log( ifelse( obs == 0, 1, obs / meanpred ) ) - (obs - meanpred), na.rm=T )
  dev.model <- 2*sum( obs *
                        log( ifelse( obs == 0, 1, obs / pred ) ) - (obs - pred) , na.rm=T )
  dev.explained <- 1 - dev.model / dev.mean
  return(c(
    dev.model,
    dev.mean,
    dev.explained) )
}

#' Internal bernoulli.dev function
#'
bernoulli.dev <- function(obs, pred){
  dev.mean <- NA
  dev.model <- NA
  dev.explained <- NA
  meanpred <- mean(obs,na.rm=T)
  dev.mean <- 2 * sum( obs*log( ifelse( meanpred == 0, 1, meanpred ) ) +
                         (1-obs)*log( ifelse( 1-meanpred == 0, 1, 1-meanpred ) ), na.rm=T )
  dev.model <- 2 * sum( obs*log( ifelse( pred == 0, 1, pred ) ) +
                          (1-obs)*log( ifelse( 1-pred == 0, 1, 1-pred ) ), na.rm=T )
  dev.explained <- 1 - dev.model / dev.mean
  return(c(
    dev.model,
    dev.mean,
    dev.explained) )
}

#' Compute PPMs
#'
#' @export
compute_ppms <- function(ppm_data_list, st_extent) {
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
  st_index <- ppm_data_list$ppm_data$date > st_extent$t.min &
              ppm_data_list$ppm_data$date < st_extent$t.max &
              ppm_data_list$ppm_data$lat > st_extent$lat.min &
              ppm_data_list$ppm_data$lat < st_extent$lat.max &
              ppm_data_list$ppm_data$lon > st_extent$lon.min &
              ppm_data_list$ppm_data$lon < st_extent$lon.max
  st_data <- ppm_data_list$ppm_data[st_index, ]

  # Add Extent of Analysis - NA & SA
  st_data$in_eoa <- FALSE
  # Add Day of Year
  st_data$DOY <- round(st_data$date * 366)

  # if lat >= 12
  ttt_index <- st_data$lat >= 12
  # by day of year
  st_data$in_eoa[ttt_index] <-
    st_data$pi.es[ttt_index] >
    ppm_data_list$eoa_es_daily$NorthAmerica[st_data$DOY[ttt_index]]

  # if lat < 12
  ttt_index <- st_data$lat < 12
  # by day of year
  st_data$in_eoa[ttt_index] <-
    st_data$pi.es[ttt_index] >
    ppm_data_list$eoa_es_daily$SouthAmerica[st_data$DOY[ttt_index]]

  # Remove data that is out of extent
  st_data <- st_data[st_data$in_eoa, ]

  # Add Binary Prediction
  # Occupied if EOA & PAT >= 0.05
  st_data$binary <- as.numeric(st_data$in_eoa & st_data$pat >= 0.05)

  # Remove rows where binary = NA (inherited from PAT = NA)
  st_data <- st_data[!is.na(st_data$binary), ]

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
    bbs <- balanced.binary.sample(binary.outcome = as.numeric(st_data$obs > 0),
                                  x = st_data$lon,
                                  y = st_data$lat,
                                  xlim = range(st_data$lon, na.rm = TRUE),
                                  ylim = range(st_data$lat, na.rm = TRUE),
                                  nx = 50,
                                  ny = 50,
                                  min.class = 0.0,
                                  neg.size = 1)

    # Index back to full vector
    sample.nindex <- c(1:nrow(st_data))[c(bbs$pos.nindex, bbs$neg.nindex)]
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
        count_stats$Spearman.abund[iii.mc] <- cor(ttt.data$pi.mu.mean,
                                                  ttt.data$obs,
                                                  method = "spearman")
        count_stats$Spearman.occ[iii.mc] <- cor(ttt.data$pi.mean,
                                                ttt.data$obs,
                                                method = "spearman")
      } # END Min Count SS
    } # END count
  } # END iii.mc

  return(list(binary_stats = binary_stats,
              occ_stats = occ_stats,
              count_stats = count_stats))
}
