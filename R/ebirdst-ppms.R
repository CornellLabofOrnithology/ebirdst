# ebirdst_ppms ----

#' eBird Status and Trends predictive performance metrics (PPMs)
#'
#' Calculate a suite of predictive performance metrics (PPMs) for the eBird
#' Status and Trends model of a given species within a spatiotemporal extent.
#'
#' @inheritParams load_raster
#' @param ext [ebirdst_extent] object (optional); the spatiotemporal extent over
#'   which to calculate the PPMs.
#' @param es_cutoff fraction between 0-1; the ensemble support cutoff to use in
#'   distinguishing zero and non-zero predictions. Optimal ensemble support
#'   cutoff values are calculated for each week during the modeling process and
#'   stored in the data package for each species. **In general, you should not
#'   specify a value for `es_cutoff` and instead allow the function to use the
#'   species-specific model-base values.**
#' @param pat_cutoff numeric between 0-1; percent above threshold.
#'
#' @details During the eBird Status and Trends modeling process, a subset of
#'   observations (the "test data") are held out from model fitting to be used
#'   for evaluating model performance. Model predictions are made for each of
#'   these observations and this function calculates a suite of predictive
#'   performance metrics (PPMs) by comparing the predictions with the observed
#'   count on the eBird checklist.
#'
#'   Three types of PPMs are calculated: binary or range-based PPMs assess the
#'   ability of model to predict range boundaries, occurrence PPMs assess the
#'   occurrence probability predictions, and abundance PPMs assess the predicted
#'   abundance. Both the occurrence and count PPMs are within-range metrics,
#'   meaning the comparison between observations and predictions is only made
#'   within the range where the species occurs.
#'
#'   Prior to calculating PPMS, the test dataset is subsampled spatiotemporally
#'   using [ebirdst_subset()]. This process is performed for 25 monte carlo
#'   iterations resulting in 25 estimates of each PPM.
#'
#' @return An `ebirdst_pppms` object containing a list of three data frames:
#'   `binary_ppms`, `occ_ppms`, and `abd_ppms`. These data frames have 25 rows
#'   corresponding to 25 Monte Carlo iterations each estimating the PPMs using a
#'   spatiotemporal subsample of the test data. Columns correspond to the
#'   different PPMS. `binary_ppms` contains binary or range-based PPMS,
#'   `occ_ppms` contains within-range occurrence probability PPMs, and
#'   `abd_ppms` contains within-range abundance PPMs. In some cases, PPMs may be
#'   missing, either because there isn't a large enough test set within the
#'   spatiotemporal extent or because average occurrence or abundance is too
#'   low. In these cases, try increasing the size of the [ebirdst_extent]
#'   object.
#'
#'   `plot()` can be called on the returned `ebirdst_ppms` object to produce a
#'   boxplot of PPMs in all three categories: Binary Occurrence, Occurrence
#'   Probability, and Abundance.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # define a spatiotemporal extent to calculate ppms over
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 42.5, ymax = 44.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' # compute predictive performance metrics
#' ppms <- ebirdst_ppms(path = path, ext = e)
#' plot(ppms)
#' }
ebirdst_ppms <- function(path, ext, es_cutoff, pat_cutoff = 1 / 7) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(is.numeric(pat_cutoff), length(pat_cutoff) == 1,
            pat_cutoff > 0, pat_cutoff < 1)
  if (!missing(ext)) {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # load configuration file
  l <- load_config(path)
  l_es_cutoff <- stats::setNames(l[["es_cutoff"]][["cutoff"]],
                                 l[["es_cutoff"]][["week"]])

  if (missing(es_cutoff)) {
    # get dynamic es cutoff
    es_cutoff <- dplyr::coalesce(l_es_cutoff, 0.75)
  } else {
    stopifnot(is.numeric(es_cutoff), length(es_cutoff) == 1,
              es_cutoff > 0, es_cutoff < 1)
    es_cutoff <- rep(es_cutoff, times = 52)
    names(es_cutoff) <- format(ebirdst::ebirdst_weeks$date, "%m-%d")
  }

  # load the test data and assign names
  preds <- load_predictions(path = path)

  # spatiotemporal subset
  if (!missing(ext)) {
    preds <- ebirdst_subset(preds, ext = ext)
  }
  if (nrow(preds) == 0) {
    warning("No predicted occurrences within spatiotemporal extent.")
    return(list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL))
  }

  # add weekly es cutoffs
  weeks <- ebirdst::ebirdst_weeks
  weeks <- weeks[format(weeks$date, "%m-%d") %in% names(es_cutoff), ]
  ws <- weeks$week_start
  we <- weeks$week_end
  if (length(es_cutoff) == 52) {
    ws[1] <- -Inf
    we[length(we)] <- Inf
  }
  preds_week <- list()
  date_frac <- preds$day_of_year / 366
  for (i in seq_along(es_cutoff)) {
    preds_week[[i]] <- preds[date_frac > ws[i] & date_frac <= we[i], ]
    preds_week[[i]][["es_cutoff"]] <- es_cutoff[i]
  }
  preds <- dplyr::bind_rows(preds_week)
  rm(preds_week, date_frac)

  # static variables
  n_mc <- 25
  # min occ sample size within range
  occ_min_ss <- 50
  occ_min_mean <- 0.01
  # min count sample size within range
  count_min_ss <- 50
  count_min_mean <- 0.25

  # define ppms
  # binary / range ppms
  binary_stat_names <- c("mc_iteration", "sample_size",
                         "mean", "auc", "pcc", "kappa",
                         "bernoulli_dev", "sensitivity", "specificity",
                         "pr_auc")
  # within range: occurrence rate ppms
  occ_stat_names <- c("mc_iteration", "sample_size",
                      "mean", "threshold", "auc", "pcc", "kappa",
                      "bernoulli_dev", "sensitivity", "specificity", "pr_auc")
  # within range: expected count ppms
  count_stat_names <- c("mc_iteration", "sample_size", "mean",
                        "poisson_dev_abd", "poisson_dev_occ",
                        "spearman_abd", "spearman_occ",
                        "spearman_count", "pearson_count_log")

  # compute monte carlo sample of ppms for spatiotemporal subset
  # split data into within range and out of range
  is_zero <- (preds$pi_es / l$FOLD_N) < preds$es_cutoff | is.na(preds$pi_es)
  zeroes <- preds[is_zero, ]
  preds <- preds[!is_zero, ]

  if (nrow(preds) == 0) {
    warning("No predicted occurrences within spatiotemporal extent.")
    return(list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL))
  }

  # false discovery rate (fdr)
  p_values <- apply(preds, 1, binom_test_p, pat_cutoff = pat_cutoff)
  p_adj <- stats::p.adjust(p_values, "fdr")
  # add binary prediction
  preds$binary <- as.numeric(p_adj < 0.01)
  # treat test data out of range with binary = 0
  if (nrow(zeroes) > 0) {
    zeroes$binary <- 0
    preds <- rbind(preds, zeroes)
  }

  # remove rows where binary is NA (inherited from pat = NA)
  # may no longer be needed
  preds <- preds[!is.na(preds$binary), ]
  if (nrow(preds) == 0) {
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
    sampled <- sample_grid(preds,
                           res = c(3000, 3000),
                           t_res = 7 / 365,
                           n = 1,
                           jitter = TRUE,
                           replace = FALSE)

    # index back to full vector
    test_sample <- preds[sampled, ]

    # binary occurrence ppms
    bs$mc_iteration[i_mc] <- i_mc
    bs$sample_size[i_mc] <- nrow(test_sample)
    bs$mean[i_mc] <- mean(as.numeric(test_sample$obs > 0))
    if (nrow(test_sample) >= occ_min_ss) {
      pa_df <- data.frame(blank = "x",
                          obs = as.numeric(test_sample$obs > 0),
                          pred = test_sample$binary)
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

      binary_curve <- precrec::evalmod(scores = test_sample$pat,
                                       labels = as.numeric(test_sample$obs > 0))
      bs$pr_auc[i_mc] <- precrec::auc(binary_curve)$aucs[2]
    }

    # within range, occurrence rate ppms
    test_inrng <- test_sample[test_sample$binary > 0, ]
    test_inrng <- test_inrng[stats::complete.cases(test_inrng$pi_median), ]

    os$mc_iteration[i_mc] <- i_mc
    os$sample_size[i_mc] <- nrow(test_inrng)
    os$mean[i_mc] <- mean(as.numeric(test_inrng$obs > 0))

    if (nrow(test_inrng) >= occ_min_ss && os$mean[i_mc] >= occ_min_mean) {
      pa_df <- data.frame(blank = "x",
                          obs = as.numeric(test_inrng$obs > 0),
                          pred = test_inrng$pi_median)
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

      occ_curve <- precrec::evalmod(scores = test_inrng$pi_median,
                                    labels = as.numeric(test_inrng$obs > 0))
      os$pr_auc[i_mc] <- precrec::auc(occ_curve)$aucs[2]
    }

    # within range, expected count ppms
    cs$mc_iteration[i_mc] <- i_mc
    cs$sample_size[i_mc] <- nrow(test_inrng)
    cs$mean[i_mc] <- mean(as.numeric(test_inrng$obs))

    if (nrow(test_inrng) >= count_min_ss && cs$mean[i_mc] >= count_min_mean) {
      # poisson deviance
      cs$poisson_dev_abd[i_mc] <- poisson_dev(obs = test_inrng$obs,
                                              pred = test_inrng$pi_mu_median)[3]

      cs$poisson_dev_occ[i_mc] <- poisson_dev(obs = test_inrng$obs,
                                              pred = test_inrng$pi_median)[3]

      # spearman's rank correlations
      cs$spearman_abd[i_mc] <- stats::cor(test_inrng$pi_mu_median,
                                          test_inrng$obs,
                                          method = "spearman")
      cs$spearman_occ[i_mc] <- stats::cor(test_inrng$pi_median,
                                          test_inrng$obs,
                                          method = "spearman")
    }

    # count ppms, based on observed occurrences
    test_cnt <- test_sample[test_sample$obs > 0, ]
    if (nrow(test_cnt) >= count_min_ss &&
        mean(test_cnt$obs, na.rm = TRUE) >= count_min_mean) {
      cs$spearman_count[i_mc] <- stats::cor(test_cnt$mu_median, test_cnt$obs,
                                            method = "spearman")
      cs$pearson_count_log[i_mc] <- stats::cor(log(test_cnt$mu_median + 1),
                                               log(test_cnt$obs + 1),
                                           method = "pearson")
    }
  }
  structure(list(binary_ppms = dplyr::tibble(type = "binary", bs),
                 occ_ppms = dplyr::tibble(type = "occurrence", os),
                 abd_ppms = dplyr::tibble(type = "abundance", cs)),
            class = "ebirdst_ppms")
}

#' @export
print.ebirdst_ppms <- function(x, ...) {
  ppm_types <- c(binary_ppms  = "Binary Range PPMs",
                 occ_ppms = "Occurrence Probability PPMs",
                 abd_ppms = "Abundance PPMs")
  stopifnot(all(names(ppm_types) %in% names(x)))

  cat("eBird Status and Trends PPMs\n")
  for (t in names(ppm_types)) {
    cat(paste0(ppm_types[t], ":\n"))
    print(x[[t]])
    cat("\n")
  }
}

#' @param x [ebirdst_ppms] object; PPMs as calculated by [ebirdst_ppms()].
#' @param ... ignored.
#'
#' @export
#' @rdname ebirdst_ppms
plot.ebirdst_ppms <- function(x, ...) {
  # transform to long
  ppms <- list()
  for (t in names(x)) {
    vars <- c("type", "auc", "pcc", "kappa", "bernoulli_dev",
              "sensitivity", "specificity", "pr_auc",
              "poisson_dev_abd", "spearman_abd",
              "spearman_count", "pearson_count_log")
    vars <- intersect(vars, names(x[[t]]))
    if (t == "binary_ppms") {
      vars <- setdiff(vars, "bernoulli_dev")
    }
    ppms[[t]] <- tidyr::pivot_longer(x[[t]][, vars],
                                     cols = -"type",
                                     names_to = "metric",
                                     values_to = "value")
  }
  ppms <- dplyr::bind_rows(ppms)
  ppms <- ppms[!is.na(ppms$value), ]

  # apply nicer labels
  ppms$label <- ppm_labels(ppms$metric)
  ppms$label <- stringr::str_replace(ppms$label, " Deviance", "\nDeviance")
  ppms$label <- stringr::str_replace(ppms$label, "Count ", "Count\n")
  ppms$label <- stringr::str_replace(ppms$label, "Abundance ", "Abundance\n")

  # construct plot for binary ppms
  ppm_b <- ppms[ppms$type == "binary", ]
  g_bin <- ggplot2::ggplot(ppm_b) +
    ggplot2::aes_string(x = "label", y = "value", group = "label") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.25) +
    ggplot2::geom_boxplot(notch = FALSE) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Binary Occurrence PPMs") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # construct plot for occurrence ppms
  ppm_o <- ppms[ppms$type == "occurrence", ]
  # negative check of bernoulli deviance
  medbde <- stats::median(ppm_o$value[ppm_o$metric == "bernoulli_dev"])
  bderep <- data.frame(type = "occurrence",
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
                  title = "Occurrence Probability PPMs") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  # if the median bernoulli deviance value is below 0, put a red x
  if (is.finite(medbde) && medbde < 0) {
    g_occ <- g_occ +
      ggplot2::geom_point(data = bderep,
                          shape = 4, size = 10, color = "red")
  }

  # construct plot for abundance ppms
  ppm_a <- ppms[ppms$type == "abundance", ]
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


# ebirdst_ppms_ts ----

#' Time series of eBird Status and Trends PPMs summarized temporally
#'
#' Calculate a time series of predictive performance metrics (PPMs) for the
#' eBird Status and Trends model. For each week or month of the year, PPMs will
#' be summarized independently to produce a time series. For further details on
#' eBird Status and Trends PPMs consult the help for [ebirdst_ppms].
#'
#' @inheritParams ebirdst_ppms
#' @param ext [ebirdst_extent] object (optional); the spatial extent over which
#'   to calculate the PPMs. Note that [ebirdst_extent] objects typically specify
#'   both a spatial and temporal extent, however, **within this function only
#'   the spatial component of the extent is used.**
#' @param summarize_by character; periods over which to summarize PPMs. PPMs can
#'   either be calculated for eBird Status and Trends weeks (as defined in
#'   [ebirdst_weeks]) or for the months of the year.
#' @param ... additional arguments passed to [ebirdst_ppms()].
#'
#' @return An `ebirdst_pppms_ts` object containing a list of three data frames:
#'   `binary_ppms`, `occ_ppms`, and `abd_ppms`. Each row of these data frames
#'   corresponds to the PPMs from one Monte Carlo iteration for a given time
#'   period. Columns correspond to the different PPMS. `binary_ppms` contains
#'   binary or range-based PPMS, `occ_ppms` contains within-range occurrence
#'   probability PPMs, and `abd_ppms` contains within-range abundance PPMs. In
#'   some cases, PPMs may be missing, either because there isn't a large enough
#'   test set within the spatiotemporal extent or because average occurrence or
#'   abundance is too low. In these cases, try increasing the size of the
#'   [ebirdst_extent] object. `plot()` can be called on the returned
#'   `ebirdst_pppms_ts` object to plot a time series of a single PPM.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # download example data
#' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' # or get the path if you already have the data downloaded
#' path <- get_species_path("example_data")
#'
#' # define a spatial extent to calculate ppms over
#' e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42.5, ymax = 44.5))
#'
#' # compute predictive performance metrics, summarized by months
#' ppms <- ebirdst_ppms_ts(path = path, ext = e, summarize_by = "months")
#'
#' # plot time series
#' # binary, kappa
#' plot(ppms, type = "binary", metric = "kappa")
#' # occurrence, sensitivity
#' plot(ppms, type = "occurrence", metric = "sensitivity")
#' #' # abundance, poisson deviance
#' plot(ppms, type = "abundance", metric = "poisson_dev_abd")
#' }
ebirdst_ppms_ts <- function(path, ext, summarize_by = c("weeks", "months"),
                            ...) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  summarize_by <- match.arg(summarize_by)

  # spatiotemporal extent
  if (missing(ext)) {
    ext <- ebirdst_extent(x = c(xmin = -180, xmax = 180,
                                ymin = -90, ymax = 90))
  } else {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }

  # set up start and end dates
  if (summarize_by == "weeks") {
    d <- ebirdst::ebirdst_weeks[, c("week_number", "date",
                                    "week_start", "week_end")]
    names(d) <- c("week_number", "date", "start", "end")
  } else if (summarize_by == "months") {
    y <- load_config(path)[["SRD_PRED_YEAR"]]
    s <- as.Date(paste(y, 1:12, "01", sep = "-"), format = "%Y-%m-%d")
    e <- s[c(2:length(s), 1)] - 1
    e <- as.Date(format(e, format = paste0(y, "-%m-%d")), format = "%Y-%m-%d")
    e <- to_srd_date(e)
    e[length(e)] <- 1
    s <- c(0, e[-length(e)])
    d <- data.frame(month = month.abb, start = s, end = e)
  }

  # calculate ppms for each period
  ppms <- list(binary_ppms = NULL, occ_ppms = NULL, abd_ppms = NULL)
  for (i in seq_len(nrow(d))) {
    # define temporal boundaries
    e <- ext
    e[["t"]][1] <- d[["start"]][i]
    e[["t"]][2] <- d[["end"]][i]

    # compute ppms
    p <- ebirdst_ppms(path = path, ext = e, ...)

    # attach date to data frames
    for (j in names(ppms)) {
      if (summarize_by == "weeks") {
        p[[j]][["week"]] <- d[["date"]][i]
        p[[j]] <- dplyr::select(p[[j]], "week", dplyr::everything())
      } else {
        p[[j]][["month"]] <- d[["month"]][i]
        p[[j]] <- dplyr::select(p[[j]], "month", dplyr::everything())
      }
      ppms[[j]] <- rbind(ppms[[j]], p[[j]])
    }
  }
  if (summarize_by == "months") {
    for (j in names(ppms)) {
      ppms[[j]][["month"]] <- factor(ppms[[j]][["month"]], month.abb)
    }
  }
  structure(ppms, class = "ebirdst_ppms_ts", summarize_by = summarize_by)
}

#' @export
print.ebirdst_ppms_ts <- function(x, ...) {
  ppm_types <- c(binary_ppms  = "Binary Range PPMs",
                 occ_ppms = "Occurrence Probability PPMs",
                 abd_ppms = "Abundance PPMs")
  stopifnot(all(names(ppm_types) %in% names(x)))

  cat("eBird Status and Trends PPMs, grouped by", attr(x, "summarize_by"), "\n")
  for (t in names(ppm_types)) {
    cat(paste0(ppm_types[t], ":\n"))
    print(x[[t]])
    cat("\n")
  }
}

#' @param x [ebirdst_ppms_ts] object; PPMs summarized by weeks or months as
#'   calculated by [ebirdst_ppms_ts()].
#' @param type character; the PPM type to plot, either a binary, occurrence, or
#'   abundance PPM can be plotted.
#' @param metric character; the specific metric to plot, the list list of
#'   possible metrics varies by PPM type:
#'   - Binary or occurrence: `auc`, `ppc`, `kappa`, `bernoulli_dev`, `sensitivity`,
#'   `specificity`
#'   - Abundance: `poisson_dev_abd`, `poisson_dev_occ`, `spearman_abd`,
#'   `spearman_occ`
#' @param ... ignored.
#'
#' @export
#' @rdname ebirdst_ppms_ts
plot.ebirdst_ppms_ts <- function(x,
                                 type = c("binary", "occurrence", "abundance"),
                                 metric = "kappa",
                                 ...) {
  type <- match.arg(type)
  stopifnot(is.character(metric), length(metric) == 1)

  # select metric to plot
  n <- dplyr::recode(type,
                     binary = "binary_ppms",
                     occurrence = "occ_ppms",
                     abundance = "abd_ppms")
  ppms <- x[[n]]
  val <- setdiff(names(ppms), c("month", "type", "mc_iteration",
                                "sample_size", "mean"))
  if (!metric %in% val) {
    stop(metric, " is not a valid ", type, " PPM. Choose one of:\n  ",
         paste(val, collapse = ", "))
  }
  if (attr(x, "summarize_by") == "months") {
    ppms <- ppms[, c("month", metric)]
  } else {
    ppms <- ppms[, c("week", metric)]
  }
  names(ppms) <- c("date", "metric")

  # plot
  if (metric == "auc") {
    metric_ylim <- c(0.5, 1)
  } else {
    metric_ylim <- c(0, 1)
  }
  # apply nicer labels
  lab <- ppm_labels(metric)
  lab <- paste0(stringr::str_to_title(type), " PPM: ", lab)

  g <- ggplot2::ggplot(ppms) +
    ggplot2::aes_string(x = "date", y = "metric", group = "date") +
    ggplot2::geom_boxplot() +
    ggplot2::ylim(metric_ylim) +
    ggplot2::labs(x = NULL, y = NULL, title = lab) +
    ggplot2::theme_light()
  if (attr(x, "summarize_by") == "weeks") {
    g <- g +
      ggplot2::scale_x_date(date_labels = "%b", date_breaks = "1 month")
  }

  suppressWarnings(print(g))
  invisible(g)
}


# internal functions ----

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

#' Binomial test for ensemble support
#'
#' @param x numeric; named numeric vector with values for `"pat"` and `"pi_es"`.
#' @param pat_cutoff numeric; percent above threshold cutoff
#'
#' @return A numeric p-value.
#'
#' @examples
#' ebirdst:::binom_test_p(c(pat = 0.1, pi_es = 75))
binom_test_p <- function(x, pat_cutoff = 1 / 10) {
  if (is.na(x["pat"]) | is.na(x["pi_es"])) {
    return(NA)
  }
  pat_pi_es <- round(as.numeric(x["pat"]) * as.numeric(x["pi_es"]), 0)
  p <- stats::binom.test(pat_pi_es, as.numeric(x["pi_es"]),
                         p = pat_cutoff,
                         alternative = "greater")
  p <- p$p.value
  return(p)
}

# apply nicer labels
ppm_labels <- function(x) {
  dplyr::case_when(
    x == "auc" ~ "ROC AUC",
    x == "pr_auc" ~ "Precision-Recall AUC",
    x == "pcc" ~ "PCC",
    x == "kappa" ~ "Kappa",
    x == "bernoulli_dev" ~ "Bernoulli Deviance",
    x == "sensitivity" ~ "Sensitivity",
    x == "specificity" ~ "Specificity",
    x == "poisson_dev_abd" ~ "Poisson Deviance",
    x == "poisson_dev_abd" ~ "Poisson Deviance",
    x == "poisson_dev_occ" ~ "Poisson Deviance",
    x == "spearman_occ" ~ "Spearman",
    x == "spearman_abd" ~ "Abundance Spearman",
    x == "spearman_count" ~ "Count Spearman",
    x == "pearson_count_log" ~ "Log Count Pearson",
    TRUE ~ x)
}
