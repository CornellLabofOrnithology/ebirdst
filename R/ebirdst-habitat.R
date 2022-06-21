#' eBird Status and Trends predictive habitat associations
#'
#' Combine the predictor importance (PI) and partial dependence (PD) data to
#' provide an estimate of the importance and directionality of the land cover
#' classes (i.e. habitat) used as covariates in the occurrence probability
#' model. **Note:** This is one of, if not the most, computationally expensive
#' operations in the package.
#'
#' @inheritParams load_pis
#' @param ext [ebirdst_extent] object; the spatiotemporal extent over which to
#'   calculate the habitat associations. Note that **temporal component of `ext`
#'   is ignored is this function**, habitat associations are always calculated
#'   for the full year.
#' @param data as an alternative to providing the `path` argument specifying the
#'   location of the data package, the data required to calculate habitat
#'   associations can be provided explicitly as a named list of three data
#'   frames: `pis` containing PI data from [load_pis()], `pds` containing PD
#'   data from [load_pds()], and `stixels` containing stixel weights in a
#'   `weight` column. All data should be provided at the stixel level,
#'   identified by the `stixel_id` column, and only those stixels appearing in
#'   the `stixels` data frame will be used. Typically stixel weights are the
#'   proportion of the focal region that the given stixel overlaps. Ignored if
#'   `path` is provided. **In most cases, users will want to avoid using these
#'   arguments and simply provide `path` instead.**
#' @param stationary_associations logical; when the habitat association should
#'   be assumed to vary throughout the year and estimates should be made for
#'   each week of the year (the default) or habitat associations should be
#'   assumed constant throughout the year and a single set of estimates made for
#'   the full year. Annual estimates should only be made when you expect the
#'   associations to be constant throughout the year, e.g. for resident species.
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
#'   - `week`: the date of the center of the week, expressed as "MM-DD". This
#'   column will be missing if `stationary_associations = TRUE`.
#'   - `importance`: the relative importance of the predictor, these values are
#'   scaled so they sum to 1 within each week.
#'   - `prob_pos_slope`: the predicted probability that the slope of the PD
#    relationship for the given predictor is positive.
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
ebirdst_habitat <- function(path, ext, data = NULL,
                            stationary_associations = FALSE) {
  stopifnot(is.logical(stationary_associations),
            length(stationary_associations) == 1)
  if (missing(path)) {
    stopifnot(is.list(data), all(c("pis", "pds", "stixels") %in% names(data)))
    pis <- data[["pis"]]
    pds <- data[["pds"]]
    stixel_coverage <- data[["stixels"]]

    col_names <- c("stixel_id", "day_of_year", "predictor", "importance")
    stopifnot(is.data.frame(pis), all(col_names %in% names(pis)))

    col_names <- c("stixel_id", "day_of_year", "predictor",
                   "predictor_value", "response")
    stopifnot(is.data.frame(pds), all(col_names %in% names(pds)))

    col_names <- c("stixel_id", "weight")
    stopifnot(is.data.frame(stixel_coverage),
              all(col_names %in% names(stixel_coverage)))
    ext <- NULL
    rm(data)
  } else {
    stopifnot(is.character(path), length(path) == 1, dir.exists(path))

    if (missing(ext)) {
      stop("A spatiotemporal extent must be provided.")
    } else {
      stopifnot(inherits(ext, "ebirdst_extent"))
      if (!sf::st_is_longlat(ext$extent)) {
        stop("Extent must provided in WGS84 lon-lat coordinates.")
      }
    }
    # remove temporal component of extent
    ext$t <- c(0, 1)
    # spatial extent polygon
    ext_poly <- ext$extent
    if (ext$type == "bbox") {
      ext_poly <- sf::st_as_sfc(ext_poly)
    }
    ext_poly <- sf::st_transform(ext_poly, crs = 4326)

    # generate stixel polygons
    if (!missing(path)) {
      stixels <- load_stixels(path = path, ext = ext)
    } else {
      stixels <- ebirdst_subset(stixels, ext = ext)
    }

    if (nrow(stixels) == 0) {
      warning("No stixels within the provided extent.")
      return(NULL)
    }

    # convert stixels to polygons
    stixels <- stixelize(stixels)
    stixels$area <- sf::st_area(stixels)
    stixels <- dplyr::select(stixels, .data$stixel_id, .data$area)

    # calculate % of stixel within focal extent
    stixels <- suppressWarnings(suppressMessages(
      sf::st_make_valid(sf::st_intersection(stixels, ext_poly))
    ))
    stixels$area_in_extent <- sf::st_area(stixels)
    stixels$weight <- as.numeric(stixels$area_in_extent / stixels$area)
    stixel_coverage <- sf::st_drop_geometry(stixels)
    stixel_coverage <- dplyr::select(stixel_coverage,
                                     .data$stixel_id,
                                     .data$weight)
    rm(stixels)

    # load pis and pds, occurrence only
    pis <- load_pis(path = path, ext = ext, model = "occurrence")
    pds <- load_pds(path = path, ext = ext, model = "occurrence")
  }

  # drop stixels not in region
  pis <- dplyr::semi_join(pis, stixel_coverage, by = "stixel_id")
  pds <- dplyr::semi_join(pds, stixel_coverage, by = "stixel_id")

  # remove intertidal if it's missing for any stixels
  it <- pis[pis$predictor == "intertidal_fs_c1_1500_pland", ]$importance
  if (any(is.na(it))) {
    pis <- pis[pis$predictor != "intertidal_fs_c1_1500_pland", ]
    pds <- pds[pds$predictor != "intertidal_fs_c1_1500_pland", ]
  }
  it <- pis[pis$predictor == "intertidal_fs_c1_1500_ed", ]$importance
  if (any(is.na(it))) {
    pis <- pis[pis$predictor != "intertidal_fs_c1_1500_ed", ]
    pds <- pds[pds$predictor != "intertidal_fs_c1_1500_ed", ]
  }
  pis <- tidyr::drop_na(pis)
  pds <- tidyr::drop_na(pds)

  if (nrow(pis) == 0 || nrow(pds) == 0) {
    warning("No stixels within the provided extent.")
    return(NULL)
  }

  # subset to cover classes
  preds <- ebirdst::ebirdst_predictors$predictor_tidy
  preds <- preds[stringr::str_detect(preds,
                                     "^(intertidal|mcd12q1|gp_rtp|astwbd|ntl)")]
  # drop variation metrics
  preds <- preds[!stringr::str_detect(preds, "_(ed|sd)$")]
  # drop unclassified
  preds <- preds[!stringr::str_detect(preds, "mcd12q1_lccs1_fs_c255")]
  pis <- pis[pis$predictor %in% preds, ]
  pds <- pds[pds$predictor %in% preds, ]

  # calculate pd slopes
  pd_slope <- pds[, c("stixel_id", "day_of_year", "predictor",
                      "predictor_value", "response")]
  rm(pds)
  pd_slope <- tidyr::nest(pd_slope, data = c("predictor_value", "response"))
  sl <- vapply(pd_slope$data, FUN = lm_slope, FUN.VALUE = numeric(1))
  # convert to binary
  pd_slope$slope <- as.numeric(sl > 0)
  pd_slope$data <- NULL
  pd_slope <- dplyr::inner_join(pd_slope, stixel_coverage, by = "stixel_id")
  rm(sl)

  # temporal smoothing of pds
  pd_smooth <- dplyr::select(pd_slope, .data$predictor,
                             x = .data$day_of_year, y = .data$slope,
                             .data$weight)
  if (isTRUE(stationary_associations)) {
    # for stationary assocations, simply calculate the weighted average
    pd_smooth <- dplyr::group_by(pd_smooth, .data$predictor)
    pd_smooth <- dplyr::summarise(pd_smooth,
                                  prob_pos_slope = weighted_mean(.data$y,
                                                                 .data$weight))
    pd_smooth <- dplyr::ungroup(pd_smooth)
  } else {
    # convert from day of year to year fraction
    pd_smooth[["x"]] <- pd_smooth[["x"]] / 366
    rm(pd_slope)
    pd_smooth <- tidyr::drop_na(pd_smooth)
    pd_smooth <- tidyr::nest(pd_smooth, data = c("x", "y", "weight"))
    pred_weeks <- ebirdst::ebirdst_weeks$week_midpoint
    pd_smooth$smooth <- lapply(pd_smooth[["data"]],
                               FUN = loess_smooth,
                               predict_to = pred_weeks,
                               na_value = NA_real_,
                               check_width = 7 / 366)
    pd_smooth$data <- NULL
    pd_smooth <- tidyr::unnest(pd_smooth, .data$smooth)
    names(pd_smooth) <- c("predictor", "week_midpoint", "prob_pos_slope")
  }

  # categorize positive and negative directionalities
  pd_smooth$prob_pos_slope <- pmin(pmax(0, pd_smooth$prob_pos_slope), 1)
  pd_smooth$direction <- NA
  pd_smooth$direction[pd_smooth$prob_pos_slope >= 0.7] <- 1
  pd_smooth$direction[pd_smooth$prob_pos_slope <= 0.3] <- -1

  # temporal smoothing of pis
  pi_smooth <- dplyr::inner_join(pis, stixel_coverage, by = "stixel_id")
  pi_smooth <- dplyr::select(pi_smooth, .data$predictor,
                             x = .data$day_of_year, y = .data$importance,
                             .data$weight)
  if (isTRUE(stationary_associations)) {
    # for stationary assocations, simply calculate the weighted average
    pi_smooth <- dplyr::group_by(pi_smooth, .data$predictor)
    pi_smooth <- dplyr::summarise(pi_smooth,
                                  importance = weighted_mean(.data$y,
                                                             .data$weight))
    pi_smooth <- dplyr::ungroup(pi_smooth)
    # scale pis
    total_importance <- sum(pi_smooth$importance, na.rm = TRUE)
    pi_smooth$importance <- pi_smooth$importance / total_importance
  } else {
    # convert from day of year to year fraction
    pi_smooth[["x"]] <- pi_smooth[["x"]] / 366
    # log transform
    pi_smooth$y <- log(pi_smooth$y + 0.001)
    pi_smooth <- tidyr::nest(pi_smooth, data = -dplyr::all_of("predictor"))
    pi_smooth$smooth <- lapply(pi_smooth[["data"]],
                               FUN = loess_smooth,
                               predict_to = pred_weeks,
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
    names(pi_smooth) <- c("predictor", "week_midpoint", "importance")
    rm(pis)
  }

  output_cols <- c("predictor", "week", "importance",
                   "prob_pos_slope", "direction")
  if (isTRUE(stationary_associations)) {
    pipd <- dplyr::inner_join(pi_smooth, pd_smooth, by = "predictor")
    output_cols <- setdiff(output_cols, "week")
  } else {
    # multiply importance and direction
    pipd <- dplyr::inner_join(pi_smooth, pd_smooth,
                              by = c("predictor", "week_midpoint"))
    # assign correct dates
    w <-  ebirdst::ebirdst_weeks
    pipd$week <- format(w$date[match(pipd$week_midpoint, w$week_midpoint)],
                        "%m-%d")
  }

  structure(pipd[, output_cols],
            class = c("ebirdst_habitat", class(pipd)),
            extent = ext)
}

#' @param x [ebirdst_habitat] object; habitat relationships as calculated by
#'   [ebirdst_habitat()].
#' @param n_habitat_types number of habitat types to include in the cake plot.
#'   The most important set of predictors will be chosen based on the maximum
#'   weekly importance value across the whole year.
#' @param ... ignored.
#'
#' @export
#' @rdname ebirdst_habitat
plot.ebirdst_habitat <- function(x, n_habitat_types = 15, ...) {
  stopifnot(is.numeric(n_habitat_types), length(n_habitat_types) == 1,
            n_habitat_types > 0)
  if (!"week" %in% names(x)) {
    stop("Habitat association plots are only possible with ",
         "stationary_associations = FALSE")
  }

  # convert to date
  x$week <- as.Date(paste(ebirdst_version()[["version_year"]],
                          x$week, sep = "-"))
  date_range <- as.Date(paste(ebirdst_version()[["version_year"]],
                              c("01-01", "12-31"), sep = "-"))

  # subset to top predictors
  x <- subset_top_predictors(x, n_habitat_types = n_habitat_types)

  # max weekly stack height
  y_max <- dplyr::group_by(x, .data$week,
                           direction = sign(.data$habitat_association))
  y_max <- dplyr::summarise(y_max,
                            height = sum(.data$habitat_association,
                                         na.rm = TRUE),
                            .groups = "drop")
  y_max <- max(abs(y_max$height), na.rm = TRUE)

  # order by importance
  x$habitat_code <- stats::reorder(x$habitat_code,
                                   x$habitat_association,
                                   FUN = function(x) {mean(abs(x))})
  # get colors
  hc <- habitat_colors[habitat_colors$habitat_code %in% x$habitat_code, ]
  hc <- dplyr::arrange(hc, .data$legend_order)

  # cake plot
  g <- ggplot2::ggplot(x) +
    ggplot2::aes_string(x = "week", y = "habitat_association",
                        fill = "habitat_code",
                        group = "habitat_code",
                        order = "habitat_code") +
    ggplot2::geom_area(colour = "white", size = 0.1) +
    ggplot2::geom_vline(xintercept = date_range[1], size = 0.25) +
    ggplot2::geom_vline(xintercept = date_range[2], size = 0.25) +
    ggplot2::geom_hline(yintercept = 0, size = 1) +
    ggplot2::scale_x_date(limits = date_range,
                          date_breaks = "1 month",
                          date_labels = "%b") +
    ggplot2::ylim(-y_max, y_max) +
    ggplot2::scale_fill_manual(values = hc$color,
                               breaks = hc$habitat_code,
                               label = hc$habitat_name) +
    ggplot2::labs(x = NULL,
                  y = "Importance (+/-)",
                  fill = "Habitat",
                  title = "Habitat associations (occurrence model)") +
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

subset_top_predictors <- function(x, n_habitat_types = 15) {
  # bring in pretty names
  lc <- dplyr::select(ebirdst::ebirdst_predictors,
                      predictor = .data$predictor_tidy,
                      .data$predictor_label,
                      .data$lc_class,
                      .data$lc_class_label)
  x <- dplyr::inner_join(lc, x, by = "predictor")

  # multiply importance and direction, fill with zeros
  x$habitat_association <- dplyr::coalesce(x$importance * x$direction, 0)

  # group roads
  x <- dplyr::group_by(x,
                       habitat_code = .data$lc_class,
                       habitat_name = .data$lc_class_label,
                       .data$week)
  x <- dplyr::summarise(x,
                        importance = sum(.data$importance, na.rm = TRUE),
                        habitat_association = sum(.data$habitat_association,
                                                  na.rm = TRUE),
                        .groups = "drop")

  # no significant associations
  if (all(is.na(x$habitat_association))) {
    warning("No significant habitat associations")
    return(x[NULL, ])
  }

  # pick top predictors based on max importance across weeks
  top_pis <- dplyr::group_by(x[abs(x$habitat_association) > 0, ],
                             .data$habitat_code)
  top_pis <- dplyr::summarise(top_pis,
                              importance = max(.data$importance),
                              .groups = "drop")
  top_pis <- dplyr::top_n(top_pis,
                          n = min(n_habitat_types, nrow(top_pis)),
                          wt = .data$importance)
  # ensure there's always 15 at most
  top_pis <- dplyr::arrange(top_pis, -.data$importance)
  top_pis <- utils::head(top_pis, 15)

  # subset to top predictors
  x[x$habitat_code %in% top_pis$habitat_code, ]
}

weighted_mean <- function(x, weight) {
  sum(x * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
}
