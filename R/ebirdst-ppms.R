#' Computes the Predictive Performance Metrics for a spatiotemporal extent
#'
#' Loads test data and ensemble support values and then calculates the
#' predictive performance metrics (PPMs) within a spatiotemporal extent defined by an
#' [ebirdst_extent] object. Use this function directly to access the computed
#' metrics, or use `plot_all_ppms()` or `plot_binary_by_time()` to summarize the
#' metrics.
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal extent to
#'   filter the data to.
#'
#' @return A list of three data frames: `binary_ppms`, `occ_ppms`, and
#'   `abd_ppms`. These data frames have 25 rows corresponding to 25 Monte Carlo
#'   iterations each estimating the PPMs using a spatiotemporal subsample of the
#'   test data. Columns correspond to the different PPMS. `binary_ppms` contains
#'   binary or range-based PPMS, `occ_ppms` contains within-range occupancy
#'   probability PPMS, and `abd_ppms` contains within-range abundance PPMs. In
#'   some cases, PPMs may be missing, either because there isn't a large enough
#'   test set within the spatiotemporal extent or because average occurrence or
#'   abundance is too low. In these cases, try increasing the size of the
#'   [ebirdst_extent] object.
#'
#' @export
#'
#' @examples
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # define a spatiotemporal extent to plot
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#' \donttest{
#' # compute predictive performance metrics
#' ppms <- compute_ppms(path = sp_path, ext = e)
#' }
compute_ppms <- function(path, ext) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # load the test data and assign names
  ppm_data <- load_test_data(path = path)
  ppm_data_raw <- load_test_data_raw(path = path)

  # add aditional rows from raw data file
  # these are assumed zeros and therefore not in test prediciton data
  ppm_data_raw <- ppm_data_raw[c("sampling_event_id", "lon", "lat", "day", "obs")]
  ppm_data_raw$date <- ppm_data_raw$day / 366
  ppm_data_raw$day <- NULL
  ppm_data_raw <- ppm_data_raw[!(ppm_data_raw$sampling_event_id %in% ppm_data$sampling_event_id), ]

  ppm_data <- dplyr::bind_rows(ppm_data, ppm_data_raw)
  rm(ppm_data_raw)

  # spatiotemporal subset
  if (!missing(ext)) {
    ppm_data <- ebirdst_subset(ppm_data, ext = ext)
  }
  if (nrow(ppm_data) == 0) {
    warning("No predicted occurrences within spatiotemporal extent.")
    return(list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL))
  }

  # static variables
  n_mc <- 25
  # min occ sample size within range
  occ_min_ss <- 50
  occ_min_mean <- 0.01
  # min count sample size within range
  count_min_ss <- 50
  count_min_mean <- 0.25
  # max count
  count_max <- NA

  # define ppms
  # binary / range ppms
  binary_stat_names <- c("mc_iteration", "sample_size",
                         "mean", "auc", "pcc", "kappa",
                         "bernoulli_dev", "sensitivity", "specificity")
  # within range: occurrence rate ppms
  occ_stat_names <- c("mc_iteration", "sample_size",
                      "mean", "threshold", "auc", "pcc", "kappa",
                      "bernoulli_dev", "sensitivity", "specificity")
  # within range: expected count ppms
  count_stat_names <- c("mc_iteration", "sample_size", "mean",
                        "poisson_dev_abd", "poisson_dev_occ",
                        "spearman_abd", "spearman_occ")

  # compute monte carlo sample of ppms for spatiotemporal subset
  # split data into within range and out of range
  ppm_data_zeroes <- ppm_data[ppm_data$pi_es < 75 | is.na(ppm_data$pi_es), ]
  ppm_data <- ppm_data[ppm_data$pi_es >= 75, ]

  if (nrow(ppm_data) == 0) {
    warning("No predicted occurrences within spatiotemporal extent.")
    return(list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL))
  }

  # false discovery rate (fdr)
  binom_test_p <- function(x) {
    if(is.na(x["pat"]) | is.na(x["pi_es"])) {
      return(NA)
    }
    pat_pi_es <- round(as.numeric(x["pat"]) * as.numeric(x["pi_es"]), 0)
    p <- stats::binom.test(pat_pi_es, as.numeric(x["pi_es"]),
                           p = (1 / 7),
                           alternative = "greater")
    p <- p$p.value
    return(p)
  }
  p_values <- apply(ppm_data, 1, binom_test_p)
  p_adj <- stats::p.adjust(p_values, "fdr")
  # add binary prediction
  ppm_data$binary <- as.numeric(p_adj < 0.01)
  # treat test data out of range with binary = 0
  if (nrow(ppm_data_zeroes) > 0) {
    ppm_data_zeroes$binary <- 0
    ppm_data <- rbind(ppm_data, ppm_data_zeroes)
  }

  # remove rows where binary is NA (inherited from pat = NA)
  # may no longer be needed
  ppm_data <- ppm_data[!is.na(ppm_data$binary), ]
  if (nrow(ppm_data) == 0) {
    warning("No predicted occurrences within spatiotemporal extent.")
    return(list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL))
  }

  # monte carlo samples for ppms: binary, spatially balanced sample
  # no oversampling of positive occurrences
  bs <- matrix(NA, n_mc, length(binary_stat_names))
  bs <- as.data.frame(bs)
  names(bs) <- binary_stat_names

  os <- matrix(NA, n_mc, length(occ_stat_names))
  os <- as.data.frame(os)
  names(os) <- occ_stat_names

  cs <- matrix(NA, n_mc, length(count_stat_names))
  cs <- as.data.frame(cs)
  names(cs) <- count_stat_names

  for (i_mc in seq_len(n_mc)) {
    # case control sampling
    sampled <- sample_case_control(ppm_data,
                                   res = c(3000, 3000),
                                   t_res = 7 / 365,
                                   n = 1,
                                   jitter = TRUE,
                                   replace = FALSE)

    # index back to full vector
    data_i <- ppm_data[sampled, ]

    # binary occupancy ppms
    bs$mc_iteration[i_mc] <- i_mc
    bs$sample_size[i_mc] <- nrow(data_i)
    bs$mean[i_mc] <- mean(as.numeric(data_i$obs > 0))
    if (nrow(data_i) >= occ_min_ss) {
      pa_df <- data.frame(blank = "x",
                          obs = as.numeric(data_i$obs > 0),
                          pred = data_i$binary)
      pa_cmx <- PresenceAbsence::cmx(pa_df, na.rm = T)
      bs$sensitivity[i_mc] <- PresenceAbsence::sensitivity(pa_cmx,
                                                           st.dev = FALSE)
      bs$specificity[i_mc] <- PresenceAbsence::specificity(pa_cmx,
                                                           st.dev = FALSE)
      bs$kappa[i_mc] <- PresenceAbsence::Kappa(pa_cmx, st.dev = FALSE)
      bs$pcc[i_mc] <- PresenceAbsence::pcc(pa_cmx, st.dev = FALSE)
      bs$auc[i_mc] <- PresenceAbsence::auc(pa_df, na.rm = TRUE,
                                           st.dev = FALSE)
      bde <- as.numeric(bernoulli_dev(obs = pa_df$obs, pred = pa_df$pred))[3]
      bs$bernoulli_dev[i_mc] <- bde
    }

    # within range, occupancy rate ppms
    data_i <- data_i[data_i$binary > 0, ]
    data_i <- data_i[stats::complete.cases(data_i$pi_mean), ]

    os$mc_iteration[i_mc] <- i_mc
    os$sample_size[i_mc] <- nrow(data_i)
    os$mean[i_mc] <- mean(as.numeric(data_i$obs > 0))

    if (nrow(data_i) >= occ_min_ss && os$mean[i_mc] >= occ_min_mean) {
      pa_df <- data.frame(blank = "x",
                          obs = as.numeric(data_i$obs > 0),
                          pred = data_i$pi_mean)
      pa_mets <- PresenceAbsence::presence.absence.accuracy(pa_df,
                                                            threshold = 0.5,
                                                            na.rm = TRUE,
                                                            st.dev = FALSE)

      # note: we use fixed estimate of the threshold, not data driven method
      #	optimal_thresh_position <- which.max(pa.metrics$Kappa)
      # this preserves the indepednece of test set
      opt_thresh_idx <- 1
      occ_thresh <- as.numeric(pa_mets$threshold[opt_thresh_idx])
      os$threshold[i_mc] <- occ_thresh
      os$pcc[i_mc] <- as.numeric(pa_mets$PCC[opt_thresh_idx])
      os$sensitivity[i_mc] <- as.numeric(pa_mets$sensitivity[opt_thresh_idx])
      os$specificity[i_mc] <- as.numeric(pa_mets$specificity[opt_thresh_idx])
      os$kappa[i_mc] <- as.numeric(pa_mets$Kappa[opt_thresh_idx])
      os$auc[i_mc] <- as.numeric(pa_mets$AUC[opt_thresh_idx])
      os$bernoulli_dev[i_mc] <- bernoulli_dev(obs = pa_df$obs,
                                              pred = pa_df$pred)[3]
    }

    # within range, expected count ppms
    # cap count at given max
    if (!is.na(count_max)) {
      data_i$obs[data_i$obs > count_max] <- count_max
    }
    cs$mc_iteration[i_mc] <- i_mc
    cs$sample_size[i_mc] <- nrow(data_i)
    cs$mean[i_mc] <- mean(as.numeric(data_i$obs))

    if(nrow(data_i) >= count_min_ss && cs$mean[i_mc] >= count_min_mean ) {
      # poisson deviance
      cs$poisson_dev_abd[i_mc] <- poisson_dev(obs = data_i$obs,
                                              pred = data_i$pi_mu_mean)[3]

      pdev <- as.numeric(poisson_dev(obs = data_i$obs,
                                     pred = data_i$pi.mean))
      cs$poisson_dev_occ[i_mc] <- poisson_dev(obs = data_i$obs,
                                              pred = data_i$pi_mean)[3]

      # spearman's rank correlations
      cs$spearman_abd[i_mc] <- stats::cor(data_i$pi_mu_mean,
                                          data_i$obs,
                                          method = "spearman")
      cs$spearman_occ[i_mc] <- stats::cor(data_i$pi_mean,
                                          data_i$obs,
                                          method = "spearman")
    }
  }
  return(list(binary_ppms = bs, occ_ppms = os, abd_ppms = cs))
}


#' Plot binary occurrence metrics by time
#'
#' For a specified number of time periods (ideally weeks or months), plots one
#' of four (Kappa, AUC, Sensitivity, Specificity) box plots. Provide an
#' `ebirdst_extent` object to see performance within a spatiotemporal extent,
#' otherwise rangewide performance will be shown.
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param metric character; the PPM to plot, eith "kappa", "auc", "sensitivity",
#'   or "specificity".
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal
#'   extent to filter the data to. The temporal component will be ignored since
#'   n_time_periods defines the temporal periods over which to calculate the
#'   PPMs.
#' @param n_time_periods integer; number of periods to divide the year into to
#'   calculate the PPMs, e.g. use 52 to divide into weeks and 12 to divide into
#'   months.
#'
#' @return Boxplot of PPM over time.
#'
#' @export
#'
#' @examples
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' \donttest{
#' # plot monthly kappa
#' plot_binary_by_time(path = sp_path, metric = "kappa", n_time_periods = 12)
#' }
plot_binary_by_time <- function(path,
                                metric = c("kappa", "auc", "sensitivity",
                                           "specificity"),
                                ext, n_time_periods = 52) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  metric <- match.arg(metric)
  stopifnot(is_integer(n_time_periods), length(n_time_periods) == 1,
            n_time_periods > 1)
  n_time_periods <- round(n_time_periods)
  if (missing(ext)) {
    ext <- ebirdst_extent(x = c(xmin = -180, xmax = 180,
                                ymin = -90, ymax = 90))
  } else {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # break into temporal units
  t_breaks <- seq(0, 1, length.out = n_time_periods + 1)
  t_dates <- t_breaks[-length(t_breaks)] + diff(t_breaks) / 2
  t_dates <- from_srd_date(t_dates)

  # calculate ppms for each
  ppms <- list(NA)
  for (i_t in seq_len(n_time_periods)) {
    e <- ebirdst_extent(x = ext$extent, t = t_breaks[c(i_t, i_t + 1)])
    ppms_i <- compute_ppms(path = path, ext = e)
    ppms_i <- ppms_i[["binary_ppms"]]
    ppms_i$date <- t_dates[[i_t]]
    ppms[[i_t]] <- ppms_i
  }
  ppms <- dplyr::bind_rows(ppms)

  stopifnot(metric %in% names(ppms))

  # plot
  if (metric == "auc") {
    metric_ylim <- c(0.5, 1)
  } else {
    metric_ylim <- c(0, 1)
  }
  metric_lab <- c(kappa = "Kappa", auc = "AUC",
                  sensitivity = "Sensitivity",
                  specificity = "Specificity")

  g <- ggplot2::ggplot(ppms) +
    ggplot2::aes_string(x = "date", y = metric, group = "date") +
    ggplot2::geom_boxplot() +
    ggplot2::ylim(metric_ylim) +
    ggplot2::labs(x = "Date", y = NULL, title = metric_lab[metric]) +
    ggplot2::scale_x_date(date_labels = "%b",
                          limits = c(as.Date("2016-01-01"),
                                     as.Date("2016-12-31")),
                          date_breaks = "1 month") +
    ggplot2::theme_light()
  suppressWarnings(print(g))
  invisible(g)
}

#' Plot all predictive performance metrics
#'
#' For a spatiotemporal extent, plots bar plots for all available predictive
#' performance metrics within three categories: Binary Occurrence, Occurrence
#' Probability, and Abundance.
#'
#' @param path character; full path to directory containing the eBird Status and
#'   Trends products for a single species.
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to filter the
#'   data to.
#'
#' @return Plot of metric box plots by category
#'
#' @export
#'
#' @examples
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -83, xmax = -82, ymin = 41, ymax = 48)
#' e <- ebirdst_extent(bb_vec, t = c("04-01", "06-30"))
#' \donttest{
#' # plot ppms within extent
#' plot_all_ppms(path = sp_path, ext = e)
#' }
plot_all_ppms <- function(path, ext) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(inherits(ext, "ebirdst_extent"))
  if (all(c(0, 1) == round(ext$t, 2))) {
    stop("Must provide temporal limits for spatiotemporal extent.")
  }

  # prepare ppms
  ppm_data <- compute_ppms(path = path, ext = ext)
  ppm_data$binary_ppms$type <- "binary"
  ppm_data$occ_ppms$type <- "occupancy"
  ppm_data$abd_ppms$type <- "abundance"
  ppm_data <- dplyr::bind_rows(ppm_data)
  ppm_data <- dplyr::select(ppm_data,
                            "type", "auc", "pcc", "kappa", "bernoulli_dev",
                            "sensitivity", "specificity",
                            "poisson_dev_abd", "spearman_abd")

  # transform to long
  ppm_data <- tidyr::gather(ppm_data, "metric", "value", -"type")
  ppm_data <- ppm_data[!(ppm_data$type == "binary" &
                           ppm_data$metric == "bernoulli_dev"), ]
  ppm_data <- ppm_data[!is.na(ppm_data$value), ]

  ppm_data$label <- dplyr::case_when(
    ppm_data$metric == "auc" ~ "AUC",
    ppm_data$metric == "pcc" ~ "PCC",
    ppm_data$metric == "kappa" ~ "Kappa",
    ppm_data$metric == "bernoulli_dev" ~ "Bernoulli\nDeviance",
    ppm_data$metric == "sensitivity" ~ "Sensitivity",
    ppm_data$metric == "specificity" ~ "Specificity",
    ppm_data$metric == "poisson_dev_abd" ~ "Poisson\nDeviance",
    ppm_data$metric == "spearman_abd" ~ "Spearman")

  # construct plot for binary ppms
  ppm_b <- ppm_data[ppm_data$type == "binary", ]
  g_bin <- ggplot2::ggplot(ppm_b) +
    ggplot2::aes_string(x = "label", y = "value", group = "label") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Binary Occupancy PPMs") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # construct plot for occupancy ppms
  ppm_o <- ppm_data[ppm_data$type == "occupancy", ]
  # negative check of bernoulli deviance
  medbde <- stats::median(ppm_o$value[ppm_o$metric == "bernoulli_dev"])
  bderep <- data.frame(type = "occupancy",
                       label = "Bernoulli\nDeviance",
                       value = 0,
                       stringsAsFactors = FALSE)
  # ggplot
  g_occ <- ggplot2::ggplot(ppm_o) +
    ggplot2::aes_string(x = "label", y = "value", group = "label") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Occupancy Probability PPMs") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  # if the median bernoulli deviance value is below 0, put a red x
  if (is.finite(medbde) && medbde < 0) {
    g_occ <- g_occ +
      ggplot2::geom_point(data = bderep,
                          shape = 4, size = 10, color = "red")
  }

  # construct plot for abundance ppms
  ppm_a <- ppm_data[ppm_data$type == "abundance", ]
  # negative check of poisson deviance
  medpde <- stats::median(ppm_a$value[ppm_a$metric == "poisson_dev_abd"])
  pderep <- data.frame(type = "abundance",
                       label = "Poisson\nDeviance",
                       value = 0,
                       stringsAsFactors = FALSE)
  # ggplot
  g_abd <- ggplot2::ggplot(ppm_a) +
    ggplot2::aes_string(x = "label", y = "value", group = "label") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Abundance PPMs") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # if the median poisson deviance value is below 0, put a red x
  if (is.finite(medpde) && medpde < 0) {
    g_abd <- g_abd +
      ggplot2::geom_point(data = pderep,
                          shape = 4, size = 10, color = "red")
  }

  three_plots <- list(g_bin, g_occ, g_abd)
  suppressWarnings(gridExtra::grid.arrange(grobs = three_plots, ncol = 3))
  invisible(three_plots)
}


#' Poisson deviance
#'
#' @param obs numeric; observed values.
#' @param pred numeric; predicted values.
#'
#' @return A named numeric vector with three elements: model deviance, mean
#'   deviance, and deviance explained.
#'
#' @examples
#' obs <- c(0, 0, 1, 3, 5, 2)
#' pred <- c(0.5, 0.1, 2.5, 3.3, 5.2, 2.5)
#' ebirdst:::poisson_dev(obs, pred)
poisson_dev <- function(obs, pred) {
  mp <- mean(obs, na.rm = TRUE)

  d_mean <- 2 * sum(obs * log(ifelse(obs == 0, 1, obs / mp)) - (obs - mp),
                    na.rm = TRUE)
  d_mod <- 2 * sum(obs * log(ifelse(obs == 0, 1, obs / pred)) - (obs - pred),
                   na.rm = TRUE)
  d_exp <- 1 - d_mod / d_mean

  c(deviance_model = d_mod, deviance_mean = d_mean, deviance_explained = d_exp)
}


#' Bernoulli deviance
#'
#' @param obs numeric; observed values.
#' @param pred numeric; predicted values.
#'
#' @return A named numeric vector with three elements: model deviance, mean
#'   deviance, and deviance explained.
#'
#' @examples
#' obs <- c(1, 1, 1, 0, 0, 0)
#' pred <- c(0.9, 0.8, 0.7, 0.3, 0.1, 0.2)
#' ebirdst:::bernoulli_dev(obs, pred)
bernoulli_dev <- function(obs, pred) {
  mp <- mean(obs, na.rm = TRUE)

  d_mean <- 2 * sum(obs * log(ifelse(mp == 0, 1, mp)) +
                      (1 - obs) * log(ifelse(1 - mp == 0, 1, 1 - mp)),
                    na.rm = TRUE)
  d_mod <- 2 * sum(obs * log(ifelse(pred == 0, 1, pred)) +
                     (1 - obs) * log(ifelse(1 - pred == 0, 1, 1 - pred)),
                   na.rm = TRUE)
  d_exp <- 1 - d_mod / d_mean

  c(deviance_model = d_mod, deviance_mean = d_mean, deviance_explained = d_exp)
}
