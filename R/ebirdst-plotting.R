#' Plot predictor importance (PI) box plots
#'
#' For a given eBird Status and Trends species, produce a box plot showing the
#' predictor importance (PI) for each of the predictors used in the occurrence
#' model. Predictors are plotted in order from highest to lowest importance.
#' Many function parameters allow for customized plots.
#'
#' @param pis data frame; predictor importance data from [load_pis()].
#' @param ext [ebirdst_extent] object; the spatiotemporal extent over which to
#'   calculate PIs. This is required, since results are less meaningful over
#'   large spatiotemporal extents.
#' @param by_cover_class logical; whether to aggregate the FRAGSTATS metrics
#'   (PLAND and ED) for the land cover classes into single values for the land
#'   cover classes.
#' @param n_top_pred integer; how many predictors to show.
#' @param pretty_names logical; whether to convert cryptic land cover codes to
#'   readable land cover class names.
#' @param plot logical; whether to plot predictor importance or just return top
#'   predictors.
#'
#' @return Plots a boxplot of predictor importance and invisibly returns a named
#'   vector of top predictors, and their median predictor importance, based on
#'   the `n_top_pred` parameter.
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
#' # load predictor importance
#' pis <- load_pis(path)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' top_pred <- plot_pis(pis, ext = e, by_cover_class = TRUE, n_top_pred = 10)
#' top_pred
#' }
plot_pis <- function(pis, ext,
                     by_cover_class = TRUE,
                     n_top_pred = 20,
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
        m[, i] <- pis[[which(lc == i)]]
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
  pis_top <- tidyr::pivot_longer(pis_top, dplyr::everything(),
                                 names_to = "predictor", values_to = "pi")

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
    y_lab <- stringr::str_glue("Relative Importance (occurrence model)")
    g <- ggplot2::ggplot(pis_top) +
      ggplot2::aes_string(x = "predictor", y = "pi") +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      ggplot2::labs(y = y_lab, x = NULL) +
      ggplot2::theme_light()
    print(g)
  }

  invisible(pi_median[top_names])
}


#' Plot partial dependency (PD) line plots
#'
#' For a given eBird Status and Trends species, produce a line plot showing the
#' partial dependence (PD) relationship for a given predictor. Two options for
#' smoothing are provided.
#'
#' @param pds data frame; partial dependence data from [load_pds()].
#' @param predictor character; single predictor name to plot PD for. For a full
#'   list of predictors, and their associated definitions, see
#'   [ebirdst_predictors].
#' @param ext [ebirdst_extent] object; the spatiotemporal extent over which to
#'   calculate PDs. This is required, since results are less meaningful over
#'   large spatiotemporal extents.
#' @param bootstrap_smooth logical; the ideal visualization of the PD data is a
#'   pointwise GAM smoothing of the individual stixel PD values. This argument
#'   specifies whether this should be done directly on the full PD dataset
#'   (`bootstrap_smooth = FALSE`) or by subsampling and bootstrapping. The
#'   latter approach deals with the randomness in the data and can be more
#'   efficient for large datasets.
#' @param show_stixel_pds logical; whether to plot the individual stixel PD
#'   values as semi-transparent lines.
#' @param show_quantiles logical; adds a band for the upper (90th) and lower
#'   (10th) quantiles of the individual stixel PD values. These are calculated
#'   using quantile regression.
#' @param n_bs int; number of GAM bootstrap iterations when estimating PD
#'   confidence intervals. Ignored if `bootstrap_smooth = FALSE`.
#' @param ss_equivalent int; when bootstrapping to estimate PD confidence
#'   intervals, this argument specifies the size of the subsample of the
#'   original data. In particular, `ss_equivalent` should be an integer
#'   representing the equivalent sampling size when averaging this number of PD
#'   estimates.
#' @param k integer; number of knots to use in the GAM when smooth the PD
#'   relationship.
#' @param ci_alpha numeric; alpha level of confidence intervals. Default is
#'   0.05.
#' @param gbm_n_trees integer; number of trees to fit in the GBM when estimating
#'   quantiles. Ignored if `show_quantiles = FALSE`. Default is 500.
#' @param ylim numeric; 2-element vector to pre-define the y-limits of plotting.
#'   In the format `c(ymin, ymax)`.
#' @param plot logical; whether to plot the PD relationships or just return
#'   data.
#'
#' @return Plots the smoothed partial dependence relationship for the specified
#'   predictor and returns a data frame of the smoothed curve with confidence
#'   intervals.
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
#' # load predictor dependence data
#' pds <- load_pds(path)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' # for testing, run with 5 bootstrap iterations for speed
#' # in practice, best to run with the default number of iterations (100)
#' pd_smooth <- plot_pds(pds, "solar_noon_diff", ext = e, n_bs = 5)
#' dplyr::glimpse(pd_smooth)
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

  # match predictor
  p <- ebirdst::ebirdst_predictors
  predictor <- stringr::str_to_lower(predictor)
  predictor <- stringr::str_replace_all(predictor, "\\.", "_")
  if (!predictor %in% p$predictor_tidy) {
    stop(paste(predictor, "is not a valid predictor variable."))
  }
  predictor_raw <- match(predictor, p$predictor_tidy)
  predictor_label <- p$predictor_label[predictor_raw]
  predictor_raw <- p$predictor_tidy[predictor_raw]

  # subset pd
  pds <- ebirdst_subset(pds, ext = ext)
  pds <- pds[pds$predictor == predictor_raw, ]
  pds <- pds[, c("stixel_id", "predictor_value", "response")]
  names(pds) <- c("stixel_id", "x", "y")

  # transform y to be relative to mean
  stix_mean <- stats::aggregate(pds["y"], list(stixel_id = pds$stixel_id),
                                mean, na.rm = TRUE)
  names(stix_mean) <- c("stixel_id", "y_mean")
  pds <- dplyr::inner_join(pds, stix_mean, by = "stixel_id")
  pds$y <- pds$y - pds$y_mean
  pds[["y_mean"]] <- NULL

  # fixed parameters
  # defines extent of pd predictions along x-axis
  # data are trimmed at both ends where they are sparser to avoid edge effects
  x_tail_level  <- 0.0
  n_preds <- 25

  # values to predict at for bootstrap iterations
  trimmed <- stats::quantile(pds$x,
                             probs = c(x_tail_level, 1 - x_tail_level),
                             na.rm = TRUE)
  nd <- data.frame(x = seq(from = trimmed[1], to = trimmed[2],
                           length = n_preds))

  # gam pointwise ci for conditional mean estimate via bootstrapping
  if (bootstrap_smooth) {
    bs_gam_pred <- matrix(NA, n_preds, n_bs)

    for (i in seq_len(n_bs)) {
      # use random sample of pd points to account for randomness in x
      # also reduces computation time
      # these ci's represent the sampling variation of an ensemble estimate
      # based on a given number of stixel pd estimates
      sample_freq <- min(ss_equivalent * n_preds / nrow(pds), 1)
      pd_sample <- dplyr::sample_frac(pds, size = sample_freq)
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
                        data = pds, gamma = 1.2)

    # calculate median and ci's
    pd_ci <- stats::predict(pd_gam, newdata = nd, se = TRUE)
    pd_ci <- data.frame(pd_median = as.numeric(pd_ci$fit),
                        pd_se = as.numeric(pd_ci$se.fit))
    pd_ci$pd_lower <- pd_ci$pd_median - 2 * pd_ci$pd_se
    pd_ci$pd_upper <- pd_ci$pd_median + 2 * pd_ci$pd_se
    pd_ci <- cbind(nd, pd_ci)
    rm(pd_gam)
  }

  # gbm quantiles
  if (isTRUE(show_quantiles)) {
    gam_ul <- gbm::gbm(y ~ x,
                       data = pds,
                       distribution = list(name = "quantile", alpha = 0.9),
                       n.trees = gbm_n_trees,
                       interaction.depth = 4,
                       shrinkage = 0.05,
                       bag.fraction = 0.5,
                       train.fraction = 1.0,
                       cv.folds = 0,
                       verbose = FALSE,
                       n.cores = 1)

    gam_ll <- gbm::gbm(y ~ x,
                       data = pds,
                       distribution = list(name = "quantile", alpha = 0.1),
                       n.trees = gbm_n_trees,
                       interaction.depth = 4,
                       shrinkage = 0.05,
                       bag.fraction = 0.5,
                       train.fraction = 1.0,
                       cv.folds = 0,
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
      g_stix <- ggplot2::geom_line(data = pds,
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
#' predictors <- c("MCD12Q1_LCCS1_FS_C1_1500_ED", "MOD44W_OIC_FS_C3_1500_ED",
#'                 "ELEV_SD")
#' ebirdst:::convert_classes(predictors, pretty = TRUE)
convert_classes <- function(x, by_cover_class = FALSE, pretty = FALSE) {
  stopifnot(is.character(x))
  stopifnot(is.logical(by_cover_class), length(by_cover_class) == 1)
  stopifnot(is.logical(pretty), length(pretty) == 1)

  x <- stringr::str_replace_all(stringr::str_to_lower(x), "\\.", "_")
  predictors_df <- ebirdst::ebirdst_predictors
  idx <- match(x, predictors_df$predictor_tidy)

  if (isTRUE(by_cover_class)) {
    if (isTRUE(pretty)) {
      y <- dplyr::coalesce(predictors_df$lc_class_label[idx], x)
    } else {
      y <- dplyr::coalesce(predictors_df$lc_class[idx],
                           predictors_df$predictor_tidy[idx],
                           x)
    }
  } else {
    if (isTRUE(pretty)) {
      y <- dplyr::coalesce(predictors_df$predictor_label[idx], x)
    } else {
      y <- dplyr::coalesce(predictors_df$predictor_tidy[idx], x)
    }
  }
  return(y)
}
