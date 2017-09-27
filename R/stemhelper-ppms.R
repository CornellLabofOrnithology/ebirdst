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
              eoa_es_daily =eoa_es_daily))
}

#' Compute PPMs
#'
#' @export
compute_ppms <- function(
  ppm.data.list,
  st.extent.list,
  n.mc = 25,
  occ.binary = T,
  occ.prob = T,
  count = T,
  # Min Occ SS within range
  occ.min.ss = 50,
  occ.min.mean = 0.01,
  # Min Occ SS within range
  count.min.ss = 50,
  count.min.mean = 0.25,
  # Set Count Max
  count.max = NA
){
  # -------------------
  # Binary / Range PPMs
  binary.stat.names <- c(
    "mc", "ss", "mean",
    "AUC", "PCC","Kappa","B.DE","Sensitivity", "Specificity")
  # Within range: Occurrence Rate PPMs
  occ.stat.names <- c(
    "mc", "ss", "mean", "threshold",
    "AUC", "PCC","Kappa","B.DE","Sensitivity", "Specificity")
  # Within range: Expected Count PPMs
  count.stat.names <- c(
    "mc","ss", "mean",
    "P.DE.abund" , "P.DE.occ", "Spearman.abund","Spearman.occ")
  # ------------------------------------------------------------------
  # Compute MC Sample of PPMs for ST Subsets
  # ------------------------------------------------------------------
  # Extract ST Subset
  st.index <-
    ppm.data.list$ppm.data$date > st.extent.list$t.min &
    ppm.data.list$ppm.data$date < st.extent.list$t.max &
    ppm.data.list$ppm.data$lat > st.extent.list$y.min &
    ppm.data.list$ppm.data$lat < st.extent.list$y.max &
    ppm.data.list$ppm.data$lon > st.extent.list$x.min &
    ppm.data.list$ppm.data$lon < st.extent.list$x.max
  st.data <- ppm.data.list$ppm.data[ st.index, ]
  # Add Extent of Analysis - NA & SA
  st.data$in.eoa <- F
  # Add Day of Year
  st.data$DOY <- round(st.data$date*366)
  # if lat >= 12
  ttt.index <- st.data$lat >= 12
  st.data$in.eoa[ttt.index] <-
    st.data$pi.es[ttt.index] >
    # by day of year
    ppm.data.list$eoa.es.daily.data$NorthAmerica[st.data$DOY[ttt.index]]
  # if lat < 12
  ttt.index <- st.data$lat < 12
  st.data$in.eoa[ttt.index] <-
    st.data$pi.es[ttt.index] >
    # by day of year
    ppm.data.list$eoa.es.daily.data$SouthAmerica[st.data$DOY[ttt.index]]
  # Remove data that is out of extent
  st.data <- st.data[ st.data$in.eoa, ]
  # Add Binary Prediction
  # Occupaied if EOA & PAT >= 0.05
  st.data$binary <- as.numeric( st.data$in.eoa & st.data$pat >= 0.05 )
  # Remove rows where binary = NA (inherited from PAT = NA)
  st.data <- st.data[ !is.na(st.data$binary), ]
  # -------------------------------------------------------------
  # Monte Carlo Samples for PPMs Binary, Spatially Balanced Sample (NO!!! oversampling +'s')
  # -------------------------------------------------------------
  binary.stats <- matrix(NA, n.mc, length(binary.stat.names))
  binary.stats <- as.data.frame(binary.stats)
  names(binary.stats) <- binary.stat.names
  occ.stats <- matrix(NA, n.mc, length(occ.stat.names))
  occ.stats <- as.data.frame(occ.stats)
  names(occ.stats) <- 	occ.stat.names
  count.stats <- matrix(NA, n.mc, length(count.stat.names))
  count.stats <- as.data.frame(count.stats)
  names(count.stats) <- count.stat.names
  for (iii.mc in 1:n.mc){
    #iii.mc <- 1
    # -------------------------------------------------------------
    # Sampling
    # -------------------------------------------------------------
    # Spatially Balanced Case Control Sampling
    # (NOTE NO!!! oversampling)
    bbs <- balanced.binary.sample(
      binary.outcome = as.numeric( st.data$obs > 0) ,
      x = st.data$lon,
      y = st.data$lat,
      xlim = range(st.data$lon, na.rm=T),
      ylim = range(st.data$lat, na.rm=T),
      nx = 50,
      ny = 50,
      min.class = 0.0,
      neg.size = 1)
    # Index back to full vector
    sample.nindex <- c(1:nrow(st.data))[ c(bbs$pos.nindex, bbs$neg.nindex) ]
    ttt.data <- st.data[sample.nindex, ]
    # Data Peek
    #bbs
    # n.cent <- plot.ppm.data.list$ppm.data(
    # 	ppm.data.list$ppm.data = ttt.data,
    # 	st.extent.list = st.extent.list,
    # 	cex = 0.5,
    # 	col = alpha("blue", 0.15) ,
    # 	pch = 19,
    # 	xlab = "Longitude",
    # 	ylab = "Latitude")
    # title(main = paste( "Selected Centroids:", n.cent) )
    # box()
    # map("state", add=T)
    # map("world", add=T)
    # -------------------------------------------------------------
    # Binary Occupancy PPMs
    # -------------------------------------------------------------
    if (occ.binary){
      binary.stats$ss[iii.mc] <- nrow(ttt.data)
      binary.stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs > 0))
      if (nrow(ttt.data) >= occ.min.ss &
          mean(as.numeric(ttt.data$obs > 0)) >= occ.min.mean ){
        # ------------
        pa.df <- data.frame(
          "space.holder.tag",
          obs= as.numeric(ttt.data$obs > 0),
          ppp= ttt.data$binary)
        pa.cmx <- cmx(pa.df, na.rm=T)
        #pa.cmx
        binary.stats$mc[iii.mc] <- iii.mc
        binary.stats$Sensitivity[iii.mc] <- sensitivity(pa.cmx, st.dev=F)
        binary.stats$Specificity[iii.mc] <- specificity(pa.cmx, st.dev=F)
        binary.stats$Kappa[iii.mc] <- Kappa(pa.cmx, st.dev=F)
        binary.stats$PCC[iii.mc] <- pcc(pa.cmx, st.dev=F)
        binary.stats$AUC[iii.mc] <- auc(pa.df, na.rm=T, st.dev=F)
        binary.stats$B.DE[iii.mc] <- as.numeric(
          bernoulli.dev(
            obs = pa.df$obs,
            pred = pa.df$ppp) )[3]
      } # END MIN SS
    } # occ.binary
    # -------------------------------------------------------------
    # Within range, Occupancy Rate PPMs
    # -------------------------------------------------------------
    if (occ.prob | count){
      # Limit ttt.data to within Range
      ttt.data <- ttt.data[ ttt.data$binary > 0, ]
    }
    if (occ.prob){
      # Remove Rows with missing pi.mean
      ttt.data <- ttt.data[ !is.na(ttt.data$pi.mean), ]
      ttt.data <- ttt.data[ !is.nan(ttt.data$pi.mean), ]
      # Min Occ SS
      occ.stats$ss[iii.mc] <- nrow(ttt.data)
      occ.stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs > 0))
      if (nrow(ttt.data) >= occ.min.ss &
          mean(as.numeric(ttt.data$obs > 0)) >= occ.min.mean ){
        # ------------
        pa2.df <- data.frame(
          "space.holder.tag",
          obs= as.numeric(ttt.data$obs > 0),
          ppp= ttt.data$pi.mean)
        pa.metrics <- presence.absence.accuracy(
          pa2.df,
          threshold = 0.5,
          #quantile(
          #pa2.df$ppp,
          #probs = seq(from=0, to=1, length=100),
          #na.rm =T),
          na.rm = T,
          st.dev = F)
        # -----------------
        # Note, we use fixed estimate of the threshold instead of the
        # data driven method,
        #	optimal_thresh_position <- which.max(pa.metrics$Kappa)
        # We do this to preserve the indepednece of test set.
        # It makes a little difference, but not much.
        # -----------------
        optimal_thresh_position <- 1
        # pa.metrics[optimal_thresh_position, ]
        # ------------------
        occ.stats$mc[iii.mc] <- iii.mc
        occ.stats$threshold[iii.mc] <-
          as.numeric(pa.metrics$threshold[optimal_thresh_position])
        occ.stats$PCC[iii.mc] <-
          as.numeric(pa.metrics$PCC[optimal_thresh_position])
        occ.stats$Sensitivity[iii.mc] <-
          as.numeric(pa.metrics$sensitivity[optimal_thresh_position])
        occ.stats$Specificity[iii.mc] <-
          as.numeric(pa.metrics$specificity[optimal_thresh_position])
        occ.stats$Kappa[iii.mc] <-
          as.numeric(pa.metrics$Kappa[optimal_thresh_position])
        occ.stats$AUC[iii.mc] <-
          as.numeric(pa.metrics$AUC[optimal_thresh_position])
        occ.stats$B.DE[iii.mc] <- as.numeric(
          bernoulli.dev(
            obs = pa2.df$obs,
            pred = pa2.df$ppp) )[3]
      } # END Min Occ SS
    } # END occ.prob
    # -------------------------------------------------------------
    # Within range, Expected Count PPMs
    # -------------------------------------------------------------
    if (count){
      # Remove Rows with missing pi.mu.mean
      ttt.data <- ttt.data[ !is.na(ttt.data$pi.mu.mean), ]
      ttt.data <- ttt.data[ !is.nan(ttt.data$pi.mu.mean), ]
      # Count Max
      if (!is.na(count.max))
        ttt.data$obs[ ttt.data$obs > count.max] <- count.max
      # Min Count SS
      count.stats$ss[iii.mc] <- nrow(ttt.data)
      count.stats$mean[iii.mc] <- mean(as.numeric(ttt.data$obs))
      # --------------
      if (nrow(ttt.data) >= count.min.ss &
          mean(as.numeric(ttt.data$obs)) >= count.min.mean ){
        # --------------
        count.stats$mc[iii.mc] <- iii.mc
        # Poisson Deviance
        pdev <- as.numeric(poisson.dev(
          obs = ttt.data$obs,
          pred = ttt.data$pi.mu.mean))
        count.stats$P.DE.abund[iii.mc] <- pdev[3]
        pdev <- as.numeric(poisson.dev(
          obs = ttt.data$obs,
          pred = ttt.data$pi.mean))
        count.stats$P.DE.occ[iii.mc] <- pdev[3]
        # Spearman's Rank Correlations
        count.stats$Spearman.abund[iii.mc] <- cor(
          ttt.data$pi.mu.mean,
          ttt.data$obs,
          method="spearman")
        count.stats$Spearman.occ[iii.mc] <- cor(
          ttt.data$pi.mean,
          ttt.data$obs,
          method="spearman")
      } # END Min Count SS
    } # END count
  } # END iii.mc
  # -------------------------------------------------------------
  return(list(
    binary.stats = binary.stats,
    occ.stats = occ.stats,
    count.stats = count.stats))
}
