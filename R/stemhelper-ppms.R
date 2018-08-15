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

  lookup.grid.cell <- function(
    xxx,
    yyy,
    ttt = NA,
    xxx.width,
    yyy.width,
    ttt.width = NA,
    jitter.cells = F ){
    # -----------------------------------
    cell.number <- rep(NA, length(xxx))
    x.ll <- y.ll <- t.ll <- NA
    nx <- ny <- nt <- NA
    if (
      length(xxx) > 0 & length(yyy) > 0 &
      length(xxx)==length(yyy) &
      xxx.width > 0 & yyy.width > 0 ){
      # Lower Left Corner
      x.ll <- min(xxx)
      y.ll <- min(yyy)
      t.ll <- min(ttt)
      # If Jittered, domain is bigger & number of grid cells increases
      if (jitter.cells){
        x.ll <- x.ll - runif(1)*xxx.width
        y.ll <- y.ll - runif(1)*yyy.width
        t.ll <- t.ll - runif(1)*ttt.width
      }
      nx <- 1 + (max(xxx) - x.ll ) %/% xxx.width
      ny <- 1 + (max(yyy) - y.ll ) %/% yyy.width
      nt <- 1 + (max(ttt) - t.ll ) %/% ttt.width
      # Assign Row, Column, Grid Cell Number
      x.number <- 1 + ( xxx - x.ll ) %/% xxx.width
      y.number <- 1 + ( yyy - y.ll ) %/% yyy.width
      t.number <- 1 + ( ttt - t.ll ) %/% ttt.width
      # Turn OFF t.number if NA
      if (all(is.na(t.number))) t.number <- 1
      cell.number <- x.number + (y.number-1)*nx + (t.number-1)*nx*ny
    }
    return( list(
      cell.number=cell.number,
      # jittered bounding box
      bb = matrix(
        c(x.ll, y.ll, t.ll,
          x.ll + nx*xxx.width, y.ll + ny*yyy.width, t.ll + nt*ttt.width),
        2, 3,
        byrow=T,
        dimnames=list(c("ll", "ur"), c("xxx", "yyy", "ttt"))),
      xxx.width = xxx.width,
      yyy.width = yyy.width,
      ttt.width = ttt.width,
      nx = nx,
      ny = ny,
      nt = nt  ))
  }

  st.grid.sampler <- function(
    xxx,
    yyy,
    ttt = NA,
    xxx.width,
    yyy.width,
    ttt.width = NA,
    jitter.cells = F,
    sample.size.per.cell = 1,
    sample.cell.probability = 1,
    replace = F ){
    # Stratified sample over Grid Cell Number
    sample_fun <- function(x, size, prob, replace){
      # Cells without samples are excluded in the tapply call -
      # if (length(x)==0) return(NA)
      # Cells with a single sample cause problems, see help(sample)
      # So, I am going to handle this situation "by hand"
      result <- rep(NA, size)
      # Flip coin to determine if cell is sampled
      if (runif(1) < prob){
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
          result <- base::sample(x=x, size=size, replace=replace)
        }
        if (length(x)>1 & replace == T ){
          result <- base::sample(x=x, size=size, replace=replace)
        }
      } # END if coin flip
      return(result)
    } # END sample_fun
    lgc <- lookup.grid.cell(
      xxx = xxx,
      yyy = yyy,
      ttt = ttt,
      xxx.width = xxx.width,
      yyy.width = yyy.width,
      ttt.width = ttt.width,
      jitter.cells = jitter.cells)
    n.index <- tapply(
      c(1:length(xxx))[!is.na(lgc$cell.number)],
      as.factor(lgc$cell.number[!is.na(lgc$cell.number)]),
      sample_fun, sample.size.per.cell, sample.cell.probability, replace)
    n.index <- plyr::rbind.fill.matrix(n.index)
    return(list(
      cell.number = lgc$cell.number,
      sample.index = n.index,
      bb = lgc$bb,
      nx = lgc$nx,
      ny = lgc$ny,
      nt = lgc$nt,
      xxx.width = lgc$xxx.width,
      yyy.width = lgc$yyy.width,
      ttt.width = lgc$ttt.width ) )
  } # END FUNCTION

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  # load configs
  e <- load_config(path)

  # load template raster
  template_raster <- raster::raster(paste(path, "/data/", e$RUN_NAME,
                                  "_srd_raster_template.tif", sep = ""))

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
  sinu <- sp::CRS(sp::proj4string(template_raster))
  st_data_prj <- sp::spTransform(st_data_sp, sinu)
  rm(st_data_sp)

  st_data_e <- raster::extract(time_stack, st_data_prj)

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

    # projected bbs
    st_data_sp <- sp::SpatialPointsDataFrame(st_data[, c("lon", "lat")],
                                             st_data,
                                             proj4string =
                                               sp::CRS("+init=epsg:4326"))
    sinu <- sp::CRS(sp::proj4string(template_raster))
    st_data_prj <- sp::spTransform(st_data_sp, sinu)

    xrange <- range(st_data_prj@coords[, 1], na.rm = TRUE)
    yrange <- range(st_data_prj@coords[, 2], na.rm = TRUE)

    xwidth <- xrange[2] - xrange[1]
    yheight <- yrange[2] - yrange[1]

    bbs <- st.grid.sampler(
      xxx = st_data_prj@coords[, 1],
      yyy = st_data_prj@coords[, 2],
      ttt = st_data$date,
      xxx.width = 10000,
      yyy.width = 10000,
      ttt.width = 7/365,
      jitter.cells = TRUE,
      sample.size.per.cell = 1,
      sample.cell.probability = 0.5,
      replace = FALSE )

    # Index back to full vector
    #sample.nindex <- c(1:nrow(st_data))[bbs$sample.index]
    ttt.data <- st_data[na.omit(bbs$sample.index), ]

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
      ttt.data <- ttt.data[ttt.data$binary > 0 & ttt.data$pi.es >= 75, ]
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
