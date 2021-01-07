#' eBird Status and Trends predictive habitat associations
#'
#' Combine the predictor importance (PI) and partial dependence (PD) data to
#' provide an estimate of the importance and directionality of the land cover
#' classes (i.e. habitat) used as model covariates. The `plot()` method can then
#' be used to produce a cake plot, a stacked area chart showing habitat
#' associations in which area indicates the importance of a given land cover
#' class and the position above or below the x-axis indicates the direction of
#' the relationship. **Note:** This is one of, if not the most, computationally
#' expensive operations in the package.
#'
#' @param path character; Full path to single species STEM results.
#' @param ext [ebirdst_extent] object; the spatiotemporal extent over which to
#'   calculate the habitat associations. Note that **temporal component of `ext`
#'   is ignored is this function**, habitat associations are always calculated
#'   for the full year.
#' @param n_predictors integer; number of land cover classes to estimate habitat
#'   associations for.
#' @param pland_only logical; For each land cover class, two FRAGSTATS metrics
#'   are used as covariates: the percent of land cover (PLAND) and the edge
#'   density (ED). By default only PLAND covariates are used to estimate habitat
#'   associates. Use `pland_only = FALSE` to use both PLAND and ED metrics.
#'
#' @return An `ebirdst_habitat` object, consisting of a data frame giving the
#'   predictor importance and directionality for each predictor for each week of
#'   the year. The columns are:
#'   - `predictor`: the name of the predictor
#'   - `date`: the week centroid expressed as a continuous value between 0-1.
#'   See [ebirdst_weeks] to convert these values to ISO dates.
#'   - `importance`: the relative importance of the predictor, these values are
#'   scaled so they sum to 1 within each week.
#'   - `direction`: the direction of the relationship, either 1 for a positive
#'   relationship, -1 for a negative relationship, or NA when the direction of
#'   the relationship is not significant.
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
#' habitat <- ebirdst_habitat(path = path, ext = e)
#' plot(habitat)
#' }
ebirdst_habitat <- function(path, ext, n_predictors = 15, pland_only = TRUE) {
  stopifnot(is.logical(by_cover_class), length(by_cover_class) == 1)
  stopifnot(is.logical(pland_only), length(pland_only) == 1)
  if (missing(ext)) {
    stop("A spatiotemporal extent must be provided.")
  } else {
    stopifnot(inherits(ext, "ebirdst_extent"))
  }
  # remove temporal component of extent
  ext$t <- c(0, 1)
  # spatial extent polygon
  ext_poly <- ext$extent
  if (ext$type == "bbox") {
    ext_poly <- sf::st_as_sfc(ext_poly)
  }
  ext_poly <- sf::st_transform(ext_poly, crs = prj_sinu)

  # generate stixel polygons
  stixels <- load_stixels(path = path)
  stixels <- ebirdst_subset(stixels, ext = ext)
  stixels <- stixelize(stixels)
  stixels <- sf::st_transform(stixels, crs = prj_sinu)
  stixels$area <- sf::st_area(stixels)
  stixels <- dplyr::select(stixels, .data$stixel_id, .data$area)

  # calculate % of stixel within focal extent
  stixels <- suppressWarnings(sf::st_intersection(stixels, ext_poly))
  stixels$area_in_extent <- sf::st_area(stixels)
  stixels$coverage <- as.numeric(stixels$area_in_extent / stixels$area)
  stixel_coverage <- sf::st_drop_geometry(stixels)
  stixel_coverage <- dplyr::select(stixel_coverage,
                                   .data$stixel_id,
                                   .data$coverage)
  rm(stixels)

  # drop any stixels that cover less than 10% of the focal extent
  stixel_coverage <- stixel_coverage[stixel_coverage$coverage >= 0.10, ]

  # load pis and pds, occurrence only
  pis <- load_pis(path = path)
  pis <- pis[pis$model == "occ", ]
  pis$model <- NULL
  pds <- load_pds(path = path, model = "occ")

  # subset
  pis <- ebirdst_subset(pis, ext = ext)
  pis <- tidyr::drop_na(pis)
  pds <- ebirdst_subset(pds, ext = ext)
  pds <- tidyr::drop_na(pds)

  # drop stixels covering less than 10% of focal area, bring in % coverage
  pis <- dplyr::inner_join(pis, stixel_coverage, by = "stixel_id")
  pds <- dplyr::inner_join(pds, stixel_coverage, by = "stixel_id")

  # pivot pis to long format
  stix_cols <- c("stixel_id", "lat", "lon", "date", "coverage")
  piv_cols <- setdiff(names(pis), stix_cols)
  pis <- tidyr::pivot_longer(pis, cols = dplyr::all_of(piv_cols),
                             names_to = "predictor",
                             values_to = "importance")

  # subset to cover classes
  preds <- ebirdst_predictors$predictor_tidy
  preds <- preds[stringr::str_detect(preds,
                                     "^(intertidal|mcd12q1|gp_rtp|astwbd|ntl)")]
  if (isTRUE(pland_only)) {
    preds <- preds[!stringr::str_detect(preds, "_ed$")]
  }
  pis <- pis[pis$predictor %in% preds, ]
  pds <- pds[pds$predictor %in% preds, ]

  # pick top predictors
  top_pis <- dplyr::group_by(pis, .data$predictor)
  top_pis <- dplyr::summarise(top_pis,
                              importance = sum(.data$importance),
                              .groups = "drop")
  top_pis <- dplyr::top_n(top_pis,
                          n = min(n_predictors, nrow(top_pis)),
                          wt = .data$importance)

  # subset to top predictors
  pis <- pis[pis$predictor %in% top_pis$predictor, ]
  pds <- pds[pds$predictor %in% top_pis$predictor, ]

  # calculate pd slopes
  pd_slope <- pds[, c("stixel_id", "date", "predictor",
                      "predictor_value", "response")]
  pd_slope <- tidyr::nest(pd_slope, data = -dplyr::all_of(c("stixel_id",
                                                            "date",
                                                            "predictor")))
  sl <- vapply(pd_slope$data, FUN = lm_slope, FUN.VALUE = numeric(1))
  # convert to binary
  pd_slope$slope <- as.numeric(sl > 0)
  pd_slope$data <- NULL
  pd_slope <- dplyr::inner_join(pd_slope, stixel_coverage, by = "stixel_id")
  rm(pds, sl)

  # aggregate across stixels with same date for speed and stability of loess
  pd_slope <- stats::aggregate(pd_slope[, c("coverage", "slope")],
                               by = pd_slope[, c("date", "predictor")],
                               FUN = mean,
                               na.rm = TRUE)

  # temporal smoothing of pds
  pd_smooth <- dplyr::select(pd_slope, .data$predictor,
                             x = .data$date, y = .data$slope,
                             weight = .data$coverage)
  pd_smooth <- tidyr::nest(pd_smooth, data = -dplyr::all_of("predictor"))
  pd_smooth$smooth <- lapply(pd_smooth[["data"]],
                             FUN = loess_smooth,
                             predict_to = ebirdst_weeks$week_midpoint,
                             na_value = NA_real_,
                             check_width = 7 / 366)
  pd_smooth$data <- NULL
  pd_smooth <- tidyr::unnest(pd_smooth, .data$smooth)
  names(pd_smooth) <- c("predictor", "date", "slope")
  rm(pd_slope)

  # categorize positive and negative directionalities
  pd_smooth$direction <- NA
  pd_smooth$direction[pd_smooth$slope >= 0.7] <- 1
  pd_smooth$direction[pd_smooth$slope <= 0.3] <- -1
  pd_smooth$slope <- NULL

  # temporal smoothing of pis
  pi_smooth <- dplyr::select(pis, .data$predictor,
                             x = .data$date, y = .data$importance,
                             weight = .data$coverage)
  # log transform
  pi_smooth$y <- log(pi_smooth$y + 0.001)
  pi_smooth <- tidyr::nest(pi_smooth, data = -dplyr::all_of("predictor"))
  pi_smooth$smooth <- lapply(pi_smooth[["data"]],
                             FUN = loess_smooth,
                             predict_to = ebirdst_weeks$week_midpoint,
                             na_value = 0,
                             check_width = 7 / 366)
  pi_smooth$data <- NULL
  pi_smooth <- tidyr::unnest(pi_smooth, .data$smooth)
  # back transform
  pi_smooth$y <- exp(pi_smooth$y)
  # scale pis
  pi_smooth <- dplyr::group_by(pi_smooth, .data$x)
  pi_smooth <- dplyr::mutate(pi_smooth,
                             y = .data$y / sum(.data$y, na.rm = TRUE))
  pi_smooth <- dplyr::ungroup(pi_smooth)
  names(pi_smooth) <- c("predictor", "date", "importance")
  rm(pis)

  # multiply importance and direction
  pipd <- dplyr::inner_join(pi_smooth, pd_smooth, by = c("predictor", "date"))

  structure(pipd, class = "ebirdst_habitat")
}

#' @param x [ebirdst_habitat] object; habitat relationships as calculated by
#'   [ebirdst_habitat()].
#' @param by_cover_class logical; whether to aggregate the FRAGSTATS metrics
#'   (PLAND and ED) for the land cover classes into single values for the land
#'   cover classes.
#' @param ... ignored.
#'
#' @export
#' @rdname ebirdst_habitat
plot.ebirdst_habitat <- function(x, by_cover_class, ....) {
  # colors for aggregated water and land cover classes
  agg_colors <- c("blue", "chartreuse", "cyan", "cyan3", "blue2", "blue4",
                  "#48D8AE", "#FCF050", "#E28373", "#C289F3", "#E6E6E6",
                  "#19AB81", "#88DA4A", "#5AB01A", "#A3B36B", "#D1BB7B",
                  "#E1D4AC", "#DCA453", "#EACA57")

  # If plotting the cover classes combined, not separately
  if(by_cover_class == TRUE) {
    # add column of class

    return_class <- function(x) {
      y <- x["predictor"]
      spls <- strsplit(y, "_")

      return(paste(spls[[1]][1], spls[[1]][2], spls[[1]][3], sep = "_"))
    }

    pipd$class <- apply(pipd, 1, return_class)

    # sum pidir by class
    pipd_agg <- stats::aggregate(pipd$pidir,
                                 by = list(pipd$class, pipd$date),
                                 FUN = sum,
                                 na.rm = TRUE)
    names(pipd_agg) <- c("predictor", "date", "pidir")

    pipd <- pipd_agg
  } else {
    agg_colors <- scales::hue_pal()(length(unique(pipd$predictor)))
  }

  # Currently, unused, but keeping it here for potential future calc
  # calculate absolute maxes of
  #absmax <- aggregate(pipd$pidir,
  #                    by = list(pipd$predictor),
  #                    FUN = function(x) { max(abs(x), na.rm = TRUE)})

  #short_set <- absmax[!is.infinite(absmax$x) & absmax$x > 0.01,]
  #pipd_short <- pipd[pipd$predictor %in% short_set$Group.1, ]

  # set final plotting object and fill with zeroes for smoothness
  pipd_short <- pipd
  pipd_short$pidir[is.na(pipd_short$pidir)] <- 0
  pipd_short$pidir[is.nan(pipd_short$pidir)] <- 0
  rm(pipd)

  # pretty the names
  ccfun <- function(x, by_cover_class) {
    convert_classes(x["predictor"],
                    by_cover_class = by_cover_class,
                    pretty = by_cover_class)
  }
  pipd_short$labels <- apply(pipd_short,
                             1,
                             FUN = ccfun,
                             by_cover_class = by_cover_class)

  if(return_data == TRUE) {
    pipd_out <- pipd_short
  }

  pipd_short$Date <- apply(pipd_short, 1, FUN = function(x) {
    strftime(as.Date(as.numeric(x["date"]) * 366,
                     origin = as.Date('2013-01-01')),
             format = "%Y-%m-%d")
  })

  pipd_short$Date <- as.Date(pipd_short$Date)
  pipd_short$date <- NULL

  # ggplot
  wave <- ggplot2::ggplot(pipd_short, ggplot2::aes(x = pipd_short$Date,
                                                   y = pipd_short$pidir,
                                                   group = pipd_short$predictor,
                                                   fill = pipd_short$predictor)) +
    ggplot2::geom_area() +
    ggplot2::geom_vline(xintercept = as.numeric(min_date)) +
    ggplot2::geom_vline(xintercept = as.numeric(max_date)) +
    ggplot2::geom_hline(yintercept = 0, size = 2) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_light() +
    ggplot2::scale_x_date(date_labels = "%b",
                          limits = c(as.Date("2013-01-01"),
                                     as.Date("2013-12-31")),
                          date_breaks = "1 month") +
    ggplot2::xlab("Date") +
    ggplot2::ylab("Association Direction (Positive/Negative)") +
    ggplot2::labs(fill = "Predictor") +
    ggplot2::theme(legend.key.size = ggplot2::unit(1, "line"))

  if(by_cover_class) {
    wave <- wave + ggplot2::scale_fill_manual(values = agg_colors,
                                              labels = pipd_short$labels,
                                              name = "Predictor")
  } else {
    wave <- wave + ggplot2::scale_fill_manual(values = agg_colors,
                                              labels = unique(pipd_short$labels),
                                              name = "Predictor")
  }

  print(wave)

  if(return_data == TRUE) {
    return(pipd_out)
  }
}

# internal functions ----

lm_slope <- function(data) {
  stopifnot(is.data.frame(data), ncol(data) == 2)
  names(data) <- c("x", "y")
  if (all(is.na(data[["y"]]))) {
    return(NA_real_)
  }
  stats::coef(stats::lm(y ~ x, data = data))[[2]]
}

loess_smooth <- function(x, predict_to, na_value = NA_real_, check_width) {
  stopifnot(is.data.frame(x), all(c("x", "y") %in% names(x)))
  stopifnot(is.numeric(predict_to))
  stopifnot(is.numeric(na_value), length(na_value) == 1)
  if (!"weight" %in% names(x)) {
    x$weight <- rep(1, times = nrow(x))
  }

  # safety check to make sure input x range is as wide as predict range
  if (min(x$x, na.rm = TRUE) > min(predict_to)) {
    if (max(x$x, na.rm = TRUE) == 1) {
      # duplicate 1 as 0
      dup <- x[x$x == 1, ]
      dup$x <- 0
      x <- rbind(x, dup)
    } else {
      # change min date to 0
      x[x$x == min(x$x, na.rm = TRUE), ]$x <- 0
    }
  }

  # need at least 2 data points to smooth
  if (sum(complete.cases(x), na.rm = TRUE) <= 1) {
    p <- return(rep(na_value, length(predict_to)))
    return(data.frame(x = predict_to, y = p))

  }
  # fit and predict
  w <- x$weight
  m <- stats::loess(formula = y ~ x,
                    degree = 1,
                    data = x,
                    weights = w)
  p <- stats::predict(m, predict_to)
  r <- data.frame(x = predict_to, y = p)

  # training data check and replace for gaps, edges, and corners
  pass <- vapply(r$x, FUN = train_check, FUN.VALUE = logical(1),
                 x_train = x$x,
                 x_range = range(predict_to),
                 check_width = check_width)
  r$y[!pass] <- na_value
  return(r)
}

train_check <- function(x, x_train, x_range, check_width) {
  # check that there is training data on either side of prediction
  # within a specified distance
  # accounting for predictions on the ends of the valid range
  if (x == x_range[1]) {
    pass <- any(x_train > x & x_train <= (x + check_width))
  } else if (x == x_range[2]) {
    pass <- any(x_train < x & x_train >= (x - check_width))
  } else {
    pass <- any(x_train > x & x_train <= (x + check_width)) &&
      any(x_train < x & x_train >= (x - check_width))
  }
  return(pass)
}
