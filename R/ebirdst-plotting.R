#' Plot predictor importances boxplot
#'
#' For all of the available predictors in a single set of species eBird
#' Status and Trends products, this function makes a bar plot of those relative
#' importances, from highest to lowest. Many function parameters allow for
#' customized plots.
#'
#' @param pis data.frame; predictor importance data from [load_pis()].
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
#'   the `n_top_pred` parameter.
#'
#' @export
#'
#' @examples
#' # download and load example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#' pis <- load_pis(sp_path)
#'
#' # define a spatiotemporal extent to plot data from
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 41.5, ymax = 43.5)
#' e <- ebirdst_extent(bb_vec, t = c("05-01", "05-31"))
#'
#' top_pred <- plot_pis(pis, ext = e, by_cover_class = TRUE, n_top_pred = 10)
#' top_pred
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

  # TODO for production, replace with ebirdst::ebirdst_predictors
  these_predictors <- ebirdst::ebirdst_predictors

  # subset
  pis <- ebirdst_subset(pis, ext = ext)
  pis <- pis[, these_predictors$predictor_tidy]

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

is_integer <- function(x) {
  is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))
}
