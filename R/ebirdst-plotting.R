#' Plot predictor importances boxplot
#'
#' For all of the available predictors in a single set of species eBird
#' Status and Trends products, this function makes a bar plot of those relative
#' importances, from highest to lowest. Many function parameters allow for
#' customized plots.
#'
#' @param pis data.frame; predictor importance data rom [load_pis()].
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to
#'   filter the data to. Required, since results are less meaningful over large
#'   spatiotemporal extents.
#' @param by_cover_class logical; whether to aggregate the four FRAGSTATS
#'   metrics for the land cover classes into single values for the land cover
#'   classes.
#' @param n_top_pred integer; how many predictors to show.
#' @param pretty_names logical; whether to convert cryptic land cover codes to
#'   readable land cover class names.
#' @param plot logical; whether to plot predictor importance or just return top
#'   predictors.
#'
#' @return Plots a boxplot of predictor importance and invisibly returns a named
#'   vector of top predictors, and their median predictor importance, based on
#'   the `n_top_pred` param.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' pis <- load_pis(sp_path)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' top_pred <- plot_pis(pis, ext = e, by_cover_class = TRUE)
#' top_pred
#'
#' }
plot_pis <- function(pis, ext,
                     by_cover_class = FALSE,
                     n_top_pred = 50,
                     pretty_names = TRUE,
                     plot = TRUE) {
  stopifnot(is.data.frame(pis))
  stopifnot(inherits(ext, "ebirdst_extent"))
  stopifnot(is.logical(by_cover_class), length(by_cover_class) == 1)
  stopifnot(is.numeric(n_top_pred), length(n_top_pred) == 1,
            n_top_pred > 1, n_top_pred <= nrow(ebirdst::ebirdst_predictors))
  stopifnot(is.logical(pretty_names), length(pretty_names) == 1)
  stopifnot(is.logical(plot), length(plot) == 1)
  if (all(c(0, 1) == round(ext$t, 2))) {
    stop("Must subset temporally, results not meaningful for full year.")
  }

  # subset
  pis <- ebirdst_subset(pis, ext = ext)
  pis <- pis[, ebirdst::ebirdst_predictors$predictor_tidy]

  # if aggregating by cover class aggregate the fragstats metrics
  if (isTRUE(by_cover_class)) {
    # find landcover classes
    lc <- convert_classes(names(pis), by_cover_class = TRUE,
                          pretty = pretty_names)
    lc_groups <- unique(lc)
    # aggregate over classes
    m <- matrix(nrow = nrow(pis), ncol = length(lc_groups))
    colnames(m) <- lc_groups
    for (i in lc_groups) {
      if (sum(lc == i) == 1) {
        m[, i] <- pis[, lc == i]
      } else {
        m[, i] <- apply(pis[, lc == i], 1, FUN = mean, na.rm = TRUE)
      }
    }
    pis <- as.data.frame(m, stringsAsFactors = FALSE)
  } else {
    names(pis) <- convert_classes(names(pis), by_cover_class = FALSE,
                                  pretty = pretty_names)
  }

  # compute median predictor importance across stixels
  pi_median <- apply(pis, 2, stats::median, na.rm = TRUE)
  pi_median <- sort(pi_median, decreasing = TRUE)

  # find the top preds based on function variable n_top_pred
  top_names <- names(pi_median)[1:min(round(n_top_pred), length(pi_median))]
  top_names <- stats::na.omit(top_names)

  # subset all values based on top_names
  pis_top <- pis[, top_names]

  # gather to long format from wide
  pis_top <- tidyr::gather(pis_top, "predictor", "pi")

  # pis have have spurious large values, NAs and NaNs
  # so clean up, trim, and check for complete cases
  pis_top$pi <- as.numeric(pis_top$pi)
  p98 <- stats::quantile(pis_top$pi, probs = 0.98, na.rm = TRUE)
  pis_top <- pis_top[pis_top$pi < p98, ]
  pis_top <- pis_top[stats::complete.cases(pis_top), ]
  pis_top$predictor <- stats::reorder(pis_top$predictor,
                                      pis_top$pi,
                                      FUN = stats::median)

  # plot
  if (isTRUE(plot)) {
    g <- ggplot2::ggplot(pis_top) +
      ggplot2::aes_string(x = "predictor", y = "pi") +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      ggplot2::labs(y = "Relative PI", x = NULL) +
      ggplot2::theme_light()
    print(g)
  }

  invisible(pi_median[top_names])
}


#' Plot partial dependency line plots
#'
#' For all of the available predictors in a single species set of eBird Status
#' and Trends products, this function makes a line plot of a single partial
#' dependency (PD), with two options for smoothing.
#'
#' @param pds data.frame; partial dependence data from `load_pds()`.
#' @param predictor character; single predictor name to plot PD for. For a full
#'   list of predictors, and their associated definitions, see
#'   [ebirdst_predictors].
#' @param ext [ebirdst_extent] object; the spatiotemporal extent to
#'   filter the data to. Required, since results are less meaningful over large
#'   spatiotemporal extents.
#' @param bootstrap_smooth logical; the ideal visualization of the PD data is a
#'   pointwise GAM smoothing of the individual stixel PD values. This argument
#'   specifies whether this should be done directly on the full PD dataset
#'   (`bootstrap_smooth = FALSE`) or by subsampling abd bootstrapping. The
#'   latter approach deals with the randomness in the data and can be more
#'   efficient for large datasets. Defaults to TRUE.
#' @param show_stixel_pds logical; whether to plot the individual stixel PD
#'   values as semi-transparent lines. Defaults to FALSE.
#' @param show_quantiles logical; adds a band for the upper (90th) and lower
#'   (10th) quantiles of the individual stixel PD values. These are calculated
#'   using quantile regression. Defaults to FALSE.
#' @param n_bs int; number of GAM bootstrap iterations when estimating PD
#'   confidence intervals. Ignored if `bootstrap_smooth = FALSE`.
#'   Default is 100.
#' @param ss_equivalent int; when bootstrapping to estimate PD confidence
#'   intervals, this argument specifies the size of the subsample of the
#'   original data. In particular, `ss_equivalent` should be an integer
#'   representing the equivalent sampling size when averaging this number of PD
#'   estimates. Defaults to 10.
#' @param k integer; number of knots to use in the GAM when smooth the PD
#'   relationship. Default is 25.
#' @param ci_alpha numeric; alpha level of confidence intervals. Default is
#'   0.05.
#' @param gbm_n_trees integer; number of trees to fit in the GBM when estimating
#'   quantiles. Ignored if `show_quantiles = FALSE`. Default is 500.
#' @param ylim numeric; 2-element vector to pre-define the y-limits of plotting.
#'   In the format c(ymin, ymax).
#' @param plot logical; whether to plot the PD relationships or just return
#'   data.
#'
#' @return Plots the smoothed partial dependence relationship and returns a data
#' frame of the smoothed curve with confidence intervals.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' pds <- load_pds(sp_path)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86.6, xmax = -82.2, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' pd_smooth <- plot_pds(pds, predictor = "time", ext = e)
#' dplyr::glimpse(pd_smooth)
#'
#' }
plot_pds <- function(pds, predictor, ext,
                     bootstrap_smooth = TRUE,
                     show_stixel_pds = FALSE,
                     show_quantiles = FALSE,
                     n_bs = 100, ss_equivalent = 10, k = 25, ci_alpha = 0.05,
                     gbm_n_trees = 500,
                     ylim = NULL,
                     plot = TRUE) {
  stopifnot(is.data.frame(pds), "predictor" %in% names(pds))
  stopifnot(is.character(predictor), length(predictor) == 1)
  stopifnot(inherits(ext, "ebirdst_extent"))
  stopifnot(is.logical(bootstrap_smooth), length(bootstrap_smooth) == 1)
  stopifnot(is.logical(show_stixel_pds), length(show_stixel_pds) == 1)
  stopifnot(is.logical(show_quantiles), length(show_quantiles) == 1)
  stopifnot(is.logical(plot), length(plot) == 1, !is.na(plot))
  stopifnot(is_integer(ss_equivalent), length(ss_equivalent) == 1,
            ss_equivalent > 0)
  stopifnot(is_integer(k), length(k) == 1, k > 0)
  stopifnot(is.numeric(ci_alpha), length(ci_alpha) == 1,
            ci_alpha > 0, ci_alpha < 0.5)
  stopifnot(is_integer(gbm_n_trees), length(gbm_n_trees) == 1, gbm_n_trees > 0)
  stopifnot(is.logical(plot), length(plot) == 1)
  if (all(c(0, 1) == round(ext$t, 2))) {
    stop("Must subset temporally, results not meaningful for full year.")
  }
  if (!is.null(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, ylim[1] < ylim[2])
  }

  # fixed parameters
  # defines extent of pd predictions along x-axis
  # data are trimmed at both ends where they are sparser to avoid edge effects
  x_tail_level  <- 0.0

  # match predictor
  predictor <- stringr::str_to_lower(predictor)
  predictor <- stringr::str_replace_all(predictor, "\\.", "_")
  if (!predictor %in% ebirdst::ebirdst_predictors$predictor_tidy) {
    stop(paste(predictor, "is not a valid predictor variable."))
  }
  predictor_raw <- match(predictor, ebirdst::ebirdst_predictors$predictor_tidy)
  predictor_raw <- ebirdst::ebirdst_predictors$predictor[predictor_raw]
  predictor_label <- convert_classes(predictor, pretty = TRUE)

  # subset pd
  pds <- ebirdst_subset(pds, ext = ext)
  pds <- pds[pds$predictor == predictor_raw, ]
  pds <- pds[!is.na(pds$y1), ]

  # transform to long format for plotting
  x_long <- pds[c("stixel_id", paste0("x", 1:50))]
  n_preds <- ncol(x_long) - 1
  x_long <- tidyr::gather(x_long, key = "key", value = "x",
                          dplyr::starts_with("x"))
  y_long <- pds[c("stixel_id", paste0("y", 1:50))]
  y_long <- tidyr::gather(y_long, key = "key", value = "y",
                          dplyr::starts_with("y"))
  pd_long <- cbind(x_long[c("stixel_id", "x")], y_long["y"])
  rm(x_long, y_long, pds)

  # transform y to be relative to mean
  stix_mean <- stats::aggregate(pd_long["y"],
                                list(stixel_id = pd_long$stixel_id),
                                mean, na.rm = TRUE)
  names(stix_mean) <- c("stixel_id", "y_mean")
  pd_long <- dplyr::inner_join(pd_long, stix_mean, by = "stixel_id")
  pd_long$y <- pd_long$y - pd_long$y_mean
  pd_long[["y_mean"]] <- NULL
  # pd_long <- dplyr::group_by(pd_long, .data$stixel_id)
  # pd_long <- dplyr::mutate(pd_long, y =  y - mean(y, na.rm = TRUE))
  # pd_long <- dplyr::ungroup(pd_long)

  # gam pointwise ci for conditional mean estimate via bootstrapping
  if (bootstrap_smooth) {
    # values to predict at for bootstrap iterations
    trimmed <- stats::quantile(pd_long$x,
                               probs = c(x_tail_level, 1 - x_tail_level),
                               na.rm = TRUE)
    nd <- data.frame(x = seq(from = trimmed[1], to = trimmed[2],
                             length = n_preds))
    bs_gam_pred <- matrix(NA, n_preds, n_bs)

    for (i in seq_len(n_bs)) {
      # use random sample of pd points to account for randomness in x
      # also reduces computation time
      # these ci's represent the sampling variation of an ensemble estimate
      # based on a given number of stixel pd estimates

      sample_freq <- ss_equivalent * n_preds / nrow(pd_long)
      pd_sample <- dplyr::sample_frac(pd_long, size = sample_freq)
      s <- mgcv::s
      pd_gam <- mgcv::gam(y ~ s(x, k = k, bs = "ds", m = 1),
                          data = pd_sample, gamma = 1.5)

      bs_gam_pred[, i] <- stats::predict(pd_gam, newdata = nd, se = FALSE)
    }

    # summarize bootstrap iterations to calc median and ci's
    pd_ci <- apply(bs_gam_pred, 1, stats::quantile,
                   probs = c(0.5, ci_alpha, 1 - ci_alpha),
                   na.rm = TRUE)
    pd_ci <- stats::setNames(as.data.frame(t(pd_ci)),
                             c("pd_median", "pd_lower", "pd_upper"))
    pd_ci <- cbind(nd, pd_ci)
    rm(pd_gam, bs_gam_pred)
  } else {
    # alternatively, use all the data
    s <- mgcv::s
    pd_gam <- mgcv::gam(y ~ s(x, k = k, bs = "ds", m = 1),
                        data = pd_long, gamma = 1.2)

    # calculate median and ci's
    pd_ci <- stats::predict(pd_gam, newdata = nd, se = TRUE)
    pd_ci <- data.frame(pd_median = pd_ci$fit, pd_se = pd_ci$se.fit)
    pd_ci$pd_lower <- pd_ci$pd_median - 2 * pd_ci$pd_se
    pd_ci$pd_upper <- pd_ci$pd_median + 2 * pd_ci$pd_se
    pd_ci <- cbind(nd, pd_ci)

    rm(pd_gam, bs_gam_pred)
  }

  # gbm quantiles
  gbm_tail_prob <- 0.1
  gbm_cv_folds <- 0
  if (isTRUE(show_quantiles)) {
    gam_ul <- gbm::gbm(y ~ x,
                       data = pd_long,
                       distribution = list(name = "quantile",
                                           alpha = (1 - gbm_tail_prob)),
                       n.trees = gbm_n_trees,
                       interaction.depth = 4,
                       shrinkage = 0.05,
                       bag.fraction = 0.5,
                       train.fraction = 1.0,
                       cv.folds = gbm_cv_folds,
                       verbose = FALSE,
                       n.cores = 1)

    gam_ll <- gbm::gbm(y ~ x,
                       data = pd_long,
                       distribution = list(name = "quantile",
                                           alpha = gbm_tail_prob),
                       n.trees = gbm_n_trees,
                       interaction.depth = 4,
                       shrinkage = 0.05,
                       bag.fraction = 0.5,
                       train.fraction = 1.0,
                       cv.folds = gbm_cv_folds,
                       verbose = FALSE,
                       n.cores = 1)

    pd_ci$pd_lower_quantile <- suppressWarnings(
      stats::predict(gam_ll, newdata = nd, n.trees = gbm_n_trees)
    )
    pd_ci$pd_upper_quantile <- suppressWarnings(
      stats::predict(gam_ul, newdata = nd, n.trees = gbm_n_trees)
    )
    rm(gam_ul, gam_ll)
  }

  if (isTRUE(plot)) {
    # y limits
    if (is.null(ylim)) {
      if (show_quantiles) {
        ylim <- c(min(c(pd_ci$pd_lower, pd_ci$pd_lower_quantile), na.rm = TRUE),
                  max(c(pd_ci$pd_upper, pd_ci$pd_upper_quantile), na.rm = TRUE))
      } else {
        ylim <- c(min(pd_ci$pd_lower, na.rm = TRUE),
                  max(pd_ci$pd_upper, na.rm = TRUE))
      }
    }

    # should stixel relationships be shown
    if (show_stixel_pds) {
      g_stix <- ggplot2::geom_line(data = pd_long,
                                   ggplot2::aes_string(x = "x", y = "y",
                                                       group = "stixel_id"),
                                   alpha = 0.2)
    } else {
      g_stix <- ggplot2::geom_blank()
    }

    # should quantiles be shown
    if (show_quantiles) {
      g_ul <- ggplot2::geom_line(
        ggplot2::aes_string(x = "x", y = "pd_upper_quantile"),
        col = "#e41a1c", linetype = "dashed")
      g_ll <- ggplot2::geom_line(
        ggplot2::aes_string(x = "x", y = "pd_lower_quantile"),
        col = "#e41a1c", linetype = "dashed")
    } else {
      g_ul <- ggplot2::geom_blank()
      g_ll <- ggplot2::geom_blank()
    }

    # main plot
    g <- ggplot2::ggplot(pd_ci) +
      ggplot2::aes_string(x = "x", y = "pd_median") +
      ggplot2::geom_hline(yintercept = 0, lwd = 0.5, col = "#000000") +
      g_stix +
      # confidence intervals
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "pd_lower",
                                               ymax = "pd_upper"),
                           fill = "#a7dBd8", alpha = 0.5) +
      # quantiles
      g_ul +
      g_ll +
      # smoothed line
      ggplot2::geom_line(col = "#fa6900", lwd = 1.5) +
      ggplot2::ylim(ylim) +
      ggplot2::labs(x = NULL, y = "Deviation E(Logit Occurrence)",
                    title = predictor_label) +
      ggplot2::theme_light()

    suppressWarnings(print(g))
  }

  invisible(pd_ci)
}


#' Converts cryptic cover class names to readable land cover names
#'
#' Internal function that converts the cryptic predictor class names to
#' readable land cover names.
#'
#' @param x character; vector of land cover variable names to convert.
#' @param by_cover_class logical; whether to replace FRAGSTATS cover class name
#'   with a name for the cover class as whole.
#' @param pretty logical; whether to convert from capital case to title case.
#'
#' @return A vector of converted names.
#'
#' @keywords internal
#'
#' @examples
#' predictors <- c("UMD_FS_C1_1500_PLAND", "MODISWATER_FS_C7_1500_LPI", "ELEV")
#' ebirdst:::convert_classes(predictors, pretty = TRUE)
convert_classes <- function(x, by_cover_class = FALSE, pretty = FALSE) {
  stopifnot(is.character(x))
  stopifnot(is.logical(by_cover_class), length(by_cover_class) == 1)
  stopifnot(is.logical(pretty), length(pretty) == 1)

  x <- stringr::str_replace_all(stringr::str_to_lower(x), "\\.", "_")
  idx <- match(x, ebirdst::ebirdst_predictors$predictor_tidy)

  if (isTRUE(by_cover_class)) {
    if (isTRUE(pretty)) {
      y <- dplyr::coalesce(ebirdst::ebirdst_predictors$lc_class_label[idx], x)
    } else {
      y <- dplyr::coalesce(ebirdst::ebirdst_predictors$lc_class[idx],
                           ebirdst::ebirdst_predictors$predictor_tidy[idx],
                           x)
    }
  } else {
    if (isTRUE(pretty)) {
      y <- dplyr::coalesce(ebirdst::ebirdst_predictors$predictor_label[idx], x)
    } else {
      y <- dplyr::coalesce(ebirdst::ebirdst_predictors$predictor_tidy[idx], x)
    }
  }
  return(y)
}

is_integer <- function(x) {
  is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))
}
