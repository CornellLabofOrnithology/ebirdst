#' eBird Status and Trends predictive habitat associations
#'
#' Combine the predictor importance (PI) and partial dependence (PD) data to
#' provide an estimate of the importance and directionality of the land cover
#' classes (i.e. habitat) used as covariates in the occurrence probability
#' model. **Note:** This is one of, if not the most, computationally expensive
#' operations in the package.
#'
#' @param path character; Full path to single species STEM results.
#' @param ext [ebirdst_extent] object; the spatiotemporal extent over which to
#'   calculate the habitat associations. Note that **temporal component of `ext`
#'   is ignored is this function**, habitat associations are always calculated
#'   for the full year.
#'
#' @details The Status and Trends models use both effort (e.g. number of
#'   observers, length of checklist) and habitat (e.g. elevation, percent forest
#'   cover) covariates; for the full list consult [ebirdst_predictors]. This
#'   function calculates habitat associations only for the following covariates
#'   that most closely represent metrics of available habitat. In all cases
#'   these are calculated within a 1.5 km radius of each checklist:
#'
#'    - Land cover: percent of each landcover class
#'    - Water cover: percent of each watercover class
#'    - Intertidal: percent cover of intertidal mudflats
#'    - Nighttime lights: total refelctance of nighttime lights
#'    - Roads: road density. There are 5 covariates distinguishing between
#'    different road types; however, these are grouped together for the sake of
#'    the habitat associations.
#'
#' The `plot()` method can be used to produce a cake plot, a stacked area chart
#' showing habitat associations in which area indicates the importance of a
#' given land cover class and the position above or below the x-axis indicates
#' the direction of the relationship.
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
#' # define a spatial extent to calculate ppms over
#' bb_vec <- c(xmin = -86, xmax = -83, ymin = 42.5, ymax = 44.5)
#' e <- ebirdst_extent(bb_vec)
#'
#' # compute habitat associations
#' habitat <- ebirdst_habitat(path = path, ext = e)
#' print(habitat)
#' # produce a cake plot
#' plot(habitat)
#' }
ebirdst_habitat <- function(path, ext) {
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))

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
  pds <- load_pds(path = path)

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
  preds <- ebirdst::ebirdst_predictors$predictor_tidy
  preds <- preds[stringr::str_detect(preds,
                                     "^(intertidal|mcd12q1|gp_rtp|astwbd|ntl)")]
  # drop variation metrics
  preds <- preds[!stringr::str_detect(preds, "_(ed|sd)$")]
  pis <- pis[pis$predictor %in% preds, ]
  pds <- pds[pds$predictor %in% preds, ]

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

  # temporal smoothing of pds
  pd_smooth <- dplyr::select(pd_slope, .data$predictor,
                             x = .data$date, y = .data$slope,
                             weight = .data$coverage)
  pd_smooth <- tidyr::drop_na(pd_smooth)
  pd_smooth <- tidyr::nest(pd_smooth, data = -dplyr::all_of("predictor"))
  pd_smooth$smooth <- lapply(pd_smooth[["data"]],
                             FUN = loess_smooth,
                             predict_to = ebirdst::ebirdst_weeks$week_midpoint,
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
                             predict_to = ebirdst::ebirdst_weeks$week_midpoint,
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

  structure(pipd,
            class = c("ebirdst_habitat", class(pipd)),
            extent = ext)
}

#' @param x [ebirdst_habitat] object; habitat relationships as calculated by
#'   [ebirdst_habitat()].
#' @param n_predictors number of predictors to include in the cake plot. The
#'   most important set of predictors will be chosen based on the maximum weekly
#'   importance value across the whole year.
#' @param date_range the range of dates for plotting; a 2-element vector of the
#'   start and end dates of the date range, provided either as dates (Date
#'   objects or strings in ISO format "YYYY-MM-DD") or numbers between 0 and 1
#'   representing the fraction of the year. When providing dates as a string,
#'   the year can be omitted (i.e. "MM-DD"). By default the full year of data
#'   are plotted.
#' @param group_roads logical; whether to aggregate the the 5 types of roads
#'   into a single class prior to plotting.
#' @param ... ignored.
#'
#' @export
#' @rdname ebirdst_habitat
plot.ebirdst_habitat <- function(x, n_predictors = 15,
                                 date_range = c(0, 1),
                                 group_roads = TRUE,
                                 ...) {
  stopifnot(is.numeric(n_predictors), length(n_predictors) == 1,
            n_predictors > 0)
  stopifnot(is.logical(group_roads), length(group_roads) == 1)

  # convert temporal extent
  date_range <- process_t_extent(date_range)
  data_date_range <- range(x$date, na.rm = TRUE)
  if (date_range[1] >= date_range[2]) {
    stop("Dates in date_range must be sequential.")
  }
  date_range[1] <- max(date_range[1], data_date_range[1])
  date_range[2] <- min(date_range[2], data_date_range[2])
  date_range <- from_srd_date(date_range)

  # bring in pretty names
  lc <- dplyr::select(ebirdst::ebirdst_predictors,
                      predictor = .data$predictor_tidy,
                      .data$predictor_label,
                      .data$lc_class,
                      .data$lc_class_label)
  x <- dplyr::inner_join(lc, x, by = "predictor")

  # multiply importance and direction, fill with zeros
  x$pi_direction <- dplyr::coalesce(x$importance * x$direction, 0)

  # group roads
  if (isTRUE(group_roads)) {
    x <- dplyr::group_by(x,
                         predictor = .data$lc_class,
                         label = .data$lc_class_label,
                         .data$date)
    x <- dplyr::summarise(x,
                          importance = sum(.data$importance),
                          pi_direction = sum(.data$pi_direction),
                          .groups = "drop")
  } else {
    x <- dplyr::select(x,
                       predictor = .data$predictor,
                       label = .data$predictor_label,
                       .data$date,
                       .data$importance,
                       .data$pi_direction)
  }

  # pick top predictors based on max importance across weeks
  top_pis <- dplyr::group_by(x, .data$predictor)
  top_pis <- dplyr::summarise(top_pis,
                              importance = max(.data$importance),
                              .groups = "drop")
  top_pis <- dplyr::top_n(top_pis,
                          n = min(n_predictors, nrow(top_pis)),
                          wt = .data$importance)

  # subset to top predictors
  x <- x[x$predictor %in% top_pis$predictor, ]

  # max weekly stack height
  y_max <- dplyr::group_by(x, .data$date,
                           direction = sign(.data$pi_direction))
  y_max <- dplyr::summarise(y_max,
                            height = sum(.data$pi_direction, na.rm = TRUE),
                            .groups = "drop")
  y_max <- max(abs(y_max$height), na.rm = TRUE)

  # assign weeks
  w <- ebirdst::ebirdst_weeks
  x$iso_date <- w$date[match(x$date, w$week_midpoint)]

  # order by importance
  x$label <- stats::reorder(x$label, x$pi_direction,
                            FUN = function(x) {mean(abs(x))})

  # cake plot
  g <- ggplot2::ggplot(x) +
    ggplot2::aes_string(x = "iso_date", y = "pi_direction",
                        group = "label", fill = "label") +
    ggplot2::geom_area(colour = "white", size = 0.1) +
    ggplot2::geom_vline(xintercept = date_range[1], size = 0.25) +
    ggplot2::geom_vline(xintercept = date_range[2], size = 0.25) +
    ggplot2::geom_hline(yintercept = 0, size = 1) +
    ggplot2::scale_x_date(limits = date_range,
                          date_breaks = "1 month",
                          date_labels = "%b") +
    ggplot2::ylim(-y_max, y_max) +
    ggplot2::scale_fill_hue() +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::labs(x = NULL,
                  y = "Importance (+/-)",
                  fill = "Predictor",
                  title = stringr::str_glue("Habitat associations ",
                                            "(occurrence model)")) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.key.size = ggplot2::unit(1, "line"))

  print(g)
  invisible(g)
}


# internal functions ----

lm_slope <- function(data) {
  stopifnot(is.data.frame(data), ncol(data) == 2)
  names(data) <- c("x", "y")
  if (all(is.na(data[["y"]])) || nrow(data) < 3) {
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

  # need at least 3 data points to smooth
  if (sum(stats::complete.cases(x), na.rm = TRUE) < 3) {
    p <- rep(na_value, length(predict_to))
    return(data.frame(x = predict_to, y = p))

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

  # fit and predict
  w <- x$weight
  suppressWarnings({
    m <- stats::loess(formula = y ~ x,
                      degree = 1,
                      data = x,
                      weights = w)
  })
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
