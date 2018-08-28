#' Plot predictor importances as barplots
#'
#' For all of the available predictors in a single species STEM result set, this
#' function makes a bar plot of those relative importances, from highest to
#' lowest. Many function parameters allow for customized plots.
#'
#' @param path character; Full path to single species STEM results.
#' @param pis data.frame; From `load_pis()`.
#' @param st_extent list; st_extent list for spatiotemporal filtering. Required,
#' as results are less meaningful over large spatiotemporal extents.
#' @param by_cover_class logical; Default is FALSE. If TRUE, aggregate Fragstats
#' for the land cover classes into single values for the land cover classes.
#' @param num_top_preds int; Integer showing how many predictors to show.
#' @param return_top logical; Default is FALSE. If TRUE, returns a vecotr of the
#' top predictors, based on the `num_top_preds` param.
#' @param pretty_names logical; Default is TRUE. Set to convert cryptic land
#' cover codes to readable land cover class names.
#' @param print_plot logical; Default is TRUE. Toggle to print plot, to allow
#' only the return the top predictors, if desired.
#'
#' @return Plots barplot and/or returns a vector of top predictors.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#' pis <- load_pis(sp_path)
#'
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' plot_pis(path = sp_path, pis = pis, st_extent = ne_extent)
#' }
plot_pis <- function(path,
                     pis,
                     st_extent,
                     by_cover_class = FALSE,
                     num_top_preds = 50,
                     return_top = FALSE,
                     pretty_names = TRUE,
                     print_plot = TRUE) {
  e <- load_config(path)

  if(num_top_preds < 2) {
    stop("num_top_preds must be greater than 1.")
  }

  if(print_plot == FALSE & return_top == FALSE) {
    stop("Both print and return params are FALSE. Nothing to do!")
  }

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  if(is.null(st_extent$t.min) | is.null(st_extent$t.max)) {
    stop("Must provide t.min and t.max as part of st_extent for this function.")
  }


  # subset for extent
  ttt_sub <- st_extent_subset(pis, st_extent)
  ttt <- ttt_sub[, e$PI_VARS]
  rm(ttt_sub)

  # if aggregating by cover class
  # aggregate the fragstats into land cover classes
  if(by_cover_class == TRUE) {
    land.cover.class.codes <- c(1:10, 12, 13, 16)
    lc.tag <- "UMD_FS_C"

    water.cover.class.codes <- c(0, 2, 3, 5, 6, 7)
    wc.tag <- "MODISWATER_FS_C"

    predictor.names <- names(ttt)

    ttt.new <- NULL
    for (iii.pred in 1:length(land.cover.class.codes)) {
      new.pred.name <- paste(lc.tag,
                             land.cover.class.codes[iii.pred],
                             "_",
                             sep = "")
      pred.nindex <- grep(new.pred.name, x = predictor.names)

      if(length(pred.nindex) == 1) {
        ttt.new <- cbind(ttt.new, ttt[, pred.nindex])
      } else {
        ttt.new <- cbind(ttt.new, apply(ttt[, pred.nindex], 1, mean, na.rm=T))
      }

      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <- paste(lc.tag,
                                             land.cover.class.codes[iii.pred],
                                             sep="")
    }
    for (iii.pred in 1:length(water.cover.class.codes)) {
      new.pred.name <- paste(wc.tag,
                             water.cover.class.codes[iii.pred],
                             "_",
                             sep = "")
      pred.nindex <- grep(new.pred.name, x = predictor.names)

      if(length(pred.nindex) == 1) {
        ttt.new <- cbind(ttt.new, ttt[, pred.nindex])
      } else {
        ttt.new <- cbind(ttt.new, apply(ttt[, pred.nindex], 1, mean, na.rm=T))
      }

      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <- paste(wc.tag,
                                             water.cover.class.codes[iii.pred],
                                             sep = "")
    }

    ttt <- ttt.new
    rm(ttt.new)
  }

  # replace names with readable
  if(pretty_names == TRUE) {
    names(ttt) <- convert_classes(names(ttt),
                                  by_cover_class = by_cover_class,
                                  pretty = by_cover_class)
  }

  # compute median
  pi_median <- apply(ttt, 2, stats::median, na.rm = T)

  # find the top preds based on function variable num_top_preds
  top_names <- names(pi_median)[order(pi_median,
                                      decreasing = T)][1:round(num_top_preds)]
  rm(pi_median)
  top_names <- stats::na.omit(top_names)

  # subset all values based on top_names
  top_pis <- ttt[, top_names]

  # munging and filtering for ggplot
  pi_stack <- utils::stack(top_pis)
  rm(top_pis)

  # PIs have have spurious large values, NAs and NaNs
  # so clean up, trim, and check for complete cases
  pi_stack$values[!is.numeric(pi_stack$values)] <- NA
  pi_stack <- pi_stack[pi_stack$values < stats::quantile(pi_stack$values,
                                                         probs = c(0.98),
                                                         na.rm = TRUE), ]
  pi_stack <- pi_stack[stats::complete.cases(pi_stack), ]

  # plot
  if(print_plot == TRUE) {
    pi_bars <- ggplot2::ggplot(pi_stack,
                               ggplot2::aes(stats::reorder(pi_stack$ind,
                                                           pi_stack$values,
                                                           FUN = stats::median),
                                            pi_stack$values)) +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      ggplot2::labs(y = "Relative PI", x = "") +
      ggplot2::theme_light()
    print(pi_bars)
  }

  if(return_top == TRUE) {
    return(top_names)
  }
}

#' Plot partial dependency as line plot
#'
#' For all of the available predictors in a single species STEM result set, this
#' function makes a line plot of a single partial dependency, with two options
#' for smoothing.
#'
#' @param pd_name character; Single predictor name from PDs (via `load_pds()`).
#' Unique predictors can be listed by calling `pds <- load_pds(path)` and
#' then calling `unique(pds$V4)`.
#' @param pds data.frame; From `load_pds()`.
#' @param st_extent list; st_extent list for spatiotemporal filtering. Required,
#' as results are less meaningful over large spatiotemporal extents.
#' @param plot_quantiles logical; Default is FALSE. Adds a band for the
#' upper (90th) and lower (10th) quantiles of the individual stixel PD values.
#' @param pointwise_pi logical; Default is TRUE. A pointwise smoothing of
#' individual stixel PD values. Ideal visualization of the PD values.
#' @param stixel_pds logical; Default is FALSE. Toggle to plot the individual
#' stixel PD values as semi-transparent lines.
#' @param k.cont.res int; Default is 25. Number of knots to use in GAM based on
#' continuining resolution of the current STEM results.
#' @param gbm.n.trees int; Default is 500. Number of trees to use in pointwise
#' GAM.
#' @param nnn.bs int; Default is 100. Daniel?
#' @param equivalent.ensemble.ss int; Default is 10. Daniel?
#' @param ci.alpha numeric; Default is 0.05. Alpha transparency of confidence
#' intervals.
#' @param mean.all.data logical; Default is FALSE. Daniel?
#' @param ylim vector pair; Opportunity to pre-define plot y-min and y-max as
#' vector pair (e.g., c(-1,1)).
#' @param print_plot logical; Default is TRUE. Set to FALSE to turn off plotting and
#' only get return of pointwise pi values.

#'
#' @return Plots barplot and returns a list containing the quantiles
#' (if plot_quantiles = TRUE) and/or the pointwise_pis upper, lower, and
#' median (if pointwise_pi = TRUE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#' pds <- load_pds(sp_path)
#'
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' plot_pds(pd_name = "TIME", pds = pds, st_extent = ne_extent)
#' }
plot_pds <- function(pd_name,
                     pds,
                     st_extent,
                     plot_quantiles = FALSE,
                     pointwise_pi = TRUE,
                     stixel_pds = FALSE,
                     k.cont.res = 25,
                     gbm.n.trees = 500,
                     nnn.bs = 100,
                     equivalent.ensemble.ss = 10,
                     ci.alpha = 0.05,
                     mean.all.data = FALSE,
                     ylim = NA,
                     print_plot = TRUE) {

  if(!(pd_name %in% unique(pds$V4))) {
    stop("Predictor name not in PDs.")
  }

  if(plot_quantiles == FALSE & pointwise_pi == FALSE & stixel_pds == FALSE) {
    stop(paste("Nothing to plot! Change one of the following to TRUE: ",
               "plot_quantiles, pointwise_pi, or stixel_pds.", sep = ""))
  }

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  if(is.null(st_extent$t.min) | is.null(st_extent$t.max)) {
    stop("Must provide t.min and t.max as part of st_extent for this function.")
  }

  # static variables
  PD_MAX_RESOLUTION <- 50
  t.ul <- NULL
  t.ll <- NULL
  t.median <- NULL
  # Defines extent of PD predictions along independent axis
  # X.tail.level of data are trimmed the both ends where data are sparser.
  # This avoids edge effects.
  x.tail.level  <- 0.0
  # Confidence level for replicate/stixel level prediction Intervals
  # Because of the relatively small sample sizes I am using 80% PIs
  gbm.tail.prob <- 0.10  # Tail probabilty [0, 0.49]
  # Number of CV folds for data driven selection of gbm's n.trees
  # parameter. Since this doesn't work on my machine I have
  # take the simpler approach of fixing the n.trees == 200
  # This seems to work well.
  gbm.cv.folds <- 0  # >= 5 computational vs statistical efficiecy
  best.iter <- gbm.n.trees
  # Number of bootstrap replicasted for Conditional Mean
  #nnn.bs <- 25
  # Number of evaluation points on x-axis / indepent var
  nd.pred.size <- PD_MAX_RESOLUTION

  return_list <- list()

  # subset based on extent
  pd_vec <- st_extent_subset(pds, st_extent)

  rm(pds)

  # Select PD Variable
  var_pd <- pd_vec[pd_vec$V4 == pd_name, ]
  rm(pd_vec)
  # Clean
  var_pd <- var_pd[!is.na(var_pd$V5), ]

  pd_name <- convert_classes(pd_name, pretty = TRUE)

  # Each Column is one replicate estimate of PD
  # 	x = x coordinate values
  # 	y = y coordinate values
  pd.x <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.y <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.mean <- rep(NA, nrow(var_pd))

  for (rid in 1:nrow(var_pd)) {
    #rid <- 100
    pd.x[, rid] <- as.numeric(
      var_pd[rid, (PD_MAX_RESOLUTION+4):(2*PD_MAX_RESOLUTION+3)])
    ttt <- as.numeric(var_pd[rid, 3:(PD_MAX_RESOLUTION+2)])
    pd.mean[rid] <- mean(ttt, na.rm=T)
    pd.y[, rid] <- ttt - pd.mean[rid]
  }

  pd.x <- as.data.frame(pd.x)
  pd.y <- as.data.frame(pd.y)
  rm(var_pd)

  # Compute Prediction Design for 1D PD
  ttt <- data.frame(x = utils::stack(pd.x)[,1],
                    y = utils::stack(pd.y)[,1])
  nd <- data.frame(x = seq(from = stats::quantile(ttt$x,
                                                  probs = x.tail.level,
                                                  na.rm = T),
                           to = stats::quantile(ttt$x,
                                                probs = 1 - x.tail.level,
                                                na.rm = T),
                           length = nd.pred.size))

  # PLOT STIXEL PD Replicates or just set up plot
  if(print_plot == TRUE) {
      if(stixel_pds) {
        if(all(is.na(ylim))) {
          ylim <- c(min(pd.y, na.rm = TRUE), max(pd.y, na.rm = TRUE))
        }

        graphics::matplot(jitter(as.matrix(pd.x), amount = 0.00),
                          pd.y,
                          ylim = ylim,
                          xlab = '',
                          ylab = "Deviation E(Logit Occurrence)",
                          type= "l",
                          lwd = 5,
                          lty = 1,
                          col = scales::alpha("black", .025))
      }
  }

  # -----------------
  # GBM Quantiles
  # -----------------
  if(plot_quantiles) {
    ttt <- data.frame(x = utils::stack(pd.x)[,1],
                      y = utils::stack(pd.y)[,1])

    ttt <- na.omit(ttt)

    d.ul <- gbm::gbm(y ~ x,
                     data = ttt,
                     distribution = list(name = "quantile",
                                         alpha = (1 - gbm.tail.prob)),
                     n.trees = gbm.n.trees,
                     interaction.depth = 4,
                     shrinkage = 0.05,
                     bag.fraction = 0.5,
                     train.fraction = 1.0,
                     cv.folds = gbm.cv.folds,
                     verbose = FALSE,
                     n.cores = 1)

    d.ll <- gbm::gbm(y ~ x,
                     data = ttt,
                     distribution = list(name = "quantile",
                                         alpha = gbm.tail.prob),
                     n.trees = gbm.n.trees,
                     interaction.depth = 4,
                     shrinkage = 0.05,
                     bag.fraction = 0.5,
                     train.fraction = 1.0,
                     cv.folds = gbm.cv.folds,
                     verbose = FALSE,
                     n.cores = 1)

    t.ul <- stats::predict(d.ul, newdata = nd, n.trees = best.iter)
    t.ll <- stats::predict(d.ll, newdata = nd, n.trees = best.iter)
    rm(d.ul, d.ll)

    poly.x <- c(nd[, 1], rev(nd[, 1]))
    poly.y <- c(t.ll, rev(t.ul))

    if(print_plot == TRUE) {

      if(stixel_pds == FALSE) {
        qymin <- min(t.ll, na.rm = TRUE)
        qymax <- max(t.ul, na.rm = TRUE)

        if(all(is.na(ylim))) {
          ylim <- c(qymin, qymax)
        }

        plot(pd.x[,1],
             pd.y[,1],
             xlab = '',
             xlim = c(min(nd, na.rm = TRUE), max(nd, na.rm = TRUE)),
             ylim = ylim,
             ylab = "Deviation E(Logit Occurrence)",
             type = "n")
      }


      graphics::polygon(poly.x,
                        poly.y,
                        col = scales::alpha("red", 0.25),
                        border = FALSE)
    }

    quantiles <- list(t.ul = t.ul, t.ll = t.ll)
    return_list$quantiles <- quantiles
  }

  # -----------------
  # GAM Pointwise CI for conditional mean estimate
  # via bootstrapping
  # -----------------
  if(pointwise_pi) {
    bs.gam.pred <- matrix(NA, nd.pred.size, nnn.bs)

    for (iii.bs in 1:nnn.bs) {
      # Take Random Sample of evaluation points
      # to account for randomness in X
      # and reduces computational time.
      # E.g. We could do this taking random
      # sample of replicates and evalution points,
      # (randomly drawing from rows and columns)
      # and set the GAM sample size at 500,
      # or the equivalent sample size when averaging
      # 10 replicate PD's each evaluated at 50 x-values.
      #
      # This suggests a nice interpretation;
      # these CI's represent the sampling variation
      # of an ensemble estimate based on a given number
      # of STIXEL PD replicates.
      #
      # equivalent.ensemble.ss <- 10
      #
      # Note that the SAME set of rows is used for every
      # column! This should be changed!
      # ------
      # Sample of STIXEL replicates
      # bs.index <- sample.int(
      # 	n = nrow(var_pd),
      # 	size = round(nrow(var_pd)*0.5),
      # 	replace = T)
      # replicate.ss <- ceiling(
      # 	equivalent.ensemble.ss*PD_MAX_RESOLUTION/length(bs.index))
      # row.index <- sample.int(
      # 	n = PD_MAX_RESOLUTION,
      # 	size = replicate.ss,
      # 	# This sets a constant fraction of available
      # 	#round(PD_MAX_RESOLUTION*0.25),
      # 	replace = F)

      rbprob <- equivalent.ensemble.ss * PD_MAX_RESOLUTION/nrow(pd.x)/ncol(pd.x)
      random.index <-  matrix((stats::rbinom(n = nrow(pd.x) * ncol(pd.x),
                                             size = 1,
                                             prob = rbprob) == 1),
                              nrow(pd.x),
                              ncol(pd.x))

      ttt <- data.frame(x = pd.x[random.index],
                        y = pd.y[random.index])

      s = mgcv::s
      d.gam <- mgcv::gam(y ~ s(x, k = k.cont.res, bs="ds", m=1),
                         data = ttt,
                         gamma = 1.5)

      bs.gam.pred[, iii.bs] <- stats::predict(d.gam, newdata = nd, se = FALSE)
      rm(d.gam)
    }

    t.ul <- apply(bs.gam.pred,
                  1,
                  stats::quantile,
                  probs = 1 - ci.alpha,
                  na.rm = TRUE)
    t.ll <- apply(bs.gam.pred,
                  1,
                  stats::quantile,
                  probs = ci.alpha,
                  na.rm = TRUE)
    t.median <- apply(bs.gam.pred,
                      1,
                      stats::quantile,
                      probs = 0.5,
                      na.rm = TRUE)

    poly.x <- c(nd[, 1], rev(nd[, 1]))
    poly.y <- c(t.ll, rev(t.ul))
    if(print_plot == TRUE) {

      if(stixel_pds == FALSE & plot_quantiles == FALSE) {
        qymin <- min(t.ll, na.rm = TRUE)
        qymax <- max(t.ul, na.rm = TRUE)

        if(all(is.na(ylim))) {
          ylim <- c(qymin, qymax)
        }

        plot(pd.x[,1],
             pd.y[,1],
             xlab = '',
             xlim = c(min(nd, na.rm = TRUE), max(nd, na.rm = TRUE)),
             ylim = ylim,
             ylab = "Deviation E(Logit Occurrence)",
             type = "n")
      }


      graphics::polygon(poly.x,
                        poly.y,
                        col = scales::alpha("blue", 0.25),
                        border = FALSE)
      graphics::lines(nd[, 1],
                      t.median,
                      col = scales::alpha("darkorange", 1.0),
                      lwd = 2 * graphics::par()$cex)
    }

    pointwise <- list(t.median = t.median, t.ul = t.ul, t.ll = t.ll)
    return_list$pointwise <- pointwise
  }

  # GAM CONDITIONAL MEAN - ALL DATA
  if(mean.all.data) {
    ttt <- data.frame(x = utils::stack(pd.x)[, 1],
                      y = utils::stack(pd.y)[, 1])

    d.gam <- mgcv::gam(y ~ s(x, k = k.cont.res, bs="ds", m=1),
                       data = ttt,
                       gamma = 1.5)

    p.gam <- stats::predict(d.gam, newdata = nd, se = TRUE)
    rm(d.gam)

    if(print_plot == TRUE) {
      graphics::polygon(x = c(nd[, 1], rev(nd[, 1])),
                        y = c(p.gam$fit + 2 * p.gam$se.fit,
                              rev(p.gam$fit - 2 * p.gam$se.fit)),
                        col = scales::alpha("lightblue", 0.25),
                        border = NA)
      graphics::lines(nd[, 1],
                      p.gam$fit,
                      col = scales::alpha("yellow", 0.75),
                      lwd = 2 * graphics::par()$cex)
    }
  }

  if(print_plot == TRUE) {
    graphics::title(pd_name, line = -2)
    graphics::abline(0, 0, col="black", lwd = 3 * graphics::par()$cex)
  }


  return(return_list)


}

#' Cake plot the PIs and PDs in a combined fasion to see directionality with
#' importance for habitat association and avoidance
#'
#' A cake plot is a stacked area chart, showing both the relative importance of
#' the land cover classes, as well as the directionality of the land cover
#' classes, by showing the stacked areas above and below a 0 line. This is one
#' of, if not the most, computationally expensive operations in the package.
#'
#' @param path character; Full path to single species STEM results.
#' @param pis data.frame; From `load_pis()`.
#' @param pds data.frame; From `load_pds()`.
#' @param st_extent list; st_extent list for spatiotemporal filtering. Required,
#' as results are less meaningful over large spatiotemporal extents.
#' @param by_cover_class logical; Default is TRUE. If TRUE, aggregate Fragstats
#' for the land cover classes into single values for the land cover classes.
#' @param pland_and_lpi_only logical; Default is TRUE. If TRUE, only the percent
#' of land cover (PLAND) and largest patch index (LPI) fragstats are used.
#' @param return_data logical; Default is FALSE. If TRUE, returns the data
#' that went into the cake plot for further manipulation.
#'
#' @return Plots a cake plot, or if return_data = TRUE, the data from the cake
#' plot as well.
#'
#' @keywords internal
#'
#' @import sp
#'
#' @examples
#' \dontrun{
#'
#' sp_path <- "path to species STEM results"
#' pis <- load_pis(sp_path)
#' pds <- load_pds(sp_path)
#'
#' ne_extent <- list(type = "rectangle",
#'                   lat.min = 40,
#'                   lat.max = 47,
#'                   lon.min = -80,
#'                   lon.max = -70,
#'                   t.min = 0.425,
#'                   t.max = 0.475)
#'
#' cake_plot(path = sp_path, pis = pis, pds = pds, st_extent = ne_extent)
#' }
cake_plot <- function(path,
                      pis,
                      pds,
                      st_extent,
                      by_cover_class = TRUE,
                      pland_and_lpi_only = TRUE,
                      return_data = FALSE) {

  ## internal loess_fit_and_predict function
  loess_fit_and_predict <- function(x, ext, input_data, type, e) {

    # internal train_check function
    train_check <- function(y, train_data, empty_val, e) {
      # static week var
      full_week <- 7/366

      # subset data to time around date
      D <- train_data
      date_sub <- D[D$date >= as.numeric(y["date"]) - full_week &
                      D$date <= as.numeric(y["date"]) + full_week, ]
      rm(D, train_data)

      if(nrow(date_sub) == 0) {
        return(NA)
      } else {
        # get max and min for comparison
        min_train <- min(date_sub$date, na.rm = TRUE)
        max_train <- max(date_sub$date, na.rm = TRUE)
        rm(date_sub)

        # end of year present problems for checking to see if there's data
        # on both sides of the prediction date, this adjusts for that
        pred_date <- as.numeric(y["date"])

        if(pred_date == e$SRD_DATE_VEC[1]) {
          pred_date <- e$SRD_DATE_VEC[1] + (full_week / 2)
        } else if(pred_date == e$SRD_DATE_VEC[length(e$SRD_DATE_VEC)]) {
          pred_date <- e$SRD_DATE_VEC[length(e$SRD_DATE_VEC)] - (full_week / 2)
        }

        if(min_train < pred_date & max_train > pred_date)  {
          # if there is data on both sides, return prediction
          return(as.numeric(y["preds"]))
        } else {
          # otherwise return empty_val
          return(empty_val)
        }
      }
    }

    # static projections
    ll <- "+init=epsg:4326"
    mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

    # for both PI and PD need to loess fit and predict for each variable
    # for uniform date set, defined here
    nd <- data.frame(date = e$SRD_DATE_VEC)

    if(type == "PI") {
      D <- input_data[, c(x, "date", "lat", "lon", "stixel_width",
                          "stixel_height")]
      names(D) <- c("predictor", "date", "lat", "lon", "stixel_width",
                    "stixel_height")
      D$predictor <- log(D$predictor + 0.001)

      empty_val <- 0
    } else {
      D <- input_data[input_data$predictor == x, ]
      D$predictor <- NULL

      empty_val <- NA
    }
    rm(input_data)

    # safety check to make sure input date range is as wide as predict dates
    if(min(D$date, na.rm = TRUE) > min(nd)) {
      if((max(D$date, na.rm = TRUE) - 1) == 0) {
        # duplicate 1 as 0
        dup <- D[D$date == 1, ]
        dup$date <- 0

        D <- rbind(D, dup)
      } else {
        # change min date to 0
        D[D$date == min(D$date, na.rm = TRUE), ]$date <- 0
      }
    }

    # make stixel polygons from centroids and widths and heights
    tdsp <- sp::SpatialPointsDataFrame(coords = D[, c("lon", "lat")],
                                       data = D,
                                       proj4string = sp::CRS(ll))

    xPlus <- tdsp$lon + (tdsp$stixel_width/2)
    yPlus <- tdsp$lat + (tdsp$stixel_height/2)
    xMinus <- tdsp$lon - (tdsp$stixel_width/2)
    yMinus <- tdsp$lat - (tdsp$stixel_height/2)

    ID <- row.names(tdsp)

    square <- cbind(xMinus, yPlus, xPlus, yPlus, xPlus,
                    yMinus, xMinus, yMinus, xMinus, yPlus)

    polys <- sp::SpatialPolygons(mapply(function(poly, id) {
      xy <- matrix(poly, ncol=2, byrow=TRUE)
      sp::Polygons(list(sp::Polygon(xy)), ID=id)
    }, split(square, row(square)), ID), proj4string=sp::CRS(ll))

    tdspolydf <- sp::SpatialPolygonsDataFrame(polys, tdsp@data)
    rm(tdsp)

    # get the full input extent
    if(ext$type == "rectangle") {
      tdsp_ext <- methods::as(raster::extent(ext$lon.min,
                                             ext$lon.max,
                                             ext$lat.min,
                                             ext$lat.max), "SpatialPolygons")
      raster::crs(tdsp_ext) <- sp::CRS(ll)
      tdsp_ext_moll <- sp::spTransform(tdsp_ext, sp::CRS(mollweide))
    } else if(ext$type == "polygon") {
      tdsp_ext_moll <- sp::spTransform(ext$polygon, sp::CRS(mollweide))
    } else {
      stop("Spatiotemporal extent type not accepted.")
    }

    # calculate area weights based on percentage of stixel covering
    # extent of interest
    tdsp_moll <- sp::spTransform(tdspolydf, sp::CRS(mollweide))
    rm(tdspolydf)
    tdsp_moll$stsqkm <- rgeos::gArea(tdsp_moll, byid = TRUE)/1000000

    rint <- raster::intersect(x = tdsp_moll, y = tdsp_ext_moll)
    rm(tdsp_moll)
    rint$intsqkm <- rgeos::gArea(rint, byid = TRUE)/1000000
    rint$refcov <- rint$intsqkm/rint$stsqkm

    # drop any stixels that cover the area less than 10%
    rint_sub <- rint[rint$refcov >= 0.10, ]
    rm(rint)

    if(type == "PI") {
      rint_d <- data.frame(predictor = rint_sub@data$predictor,
                           date = rint_sub@data$date)

      y_name <- "predictor"
    } else {
      rint_d <- data.frame(slope = rint_sub@data$slope,
                           date = rint_sub@data$date)

      y_name <- "slope"
    }

    if(!all(is.nan(rint_d[, c(y_name)])) &
       sum(!is.nan(rint_d[, c(y_name)]), na.rm = TRUE) > 1) {
      d.loess <- stats::loess(formula = stats::as.formula(paste(y_name,
                                                                " ~ date",
                                                                sep = "")),
                              degree = 1,
                              data = rint_d,
                              weights = rint_sub@data$refcov)
      rm(rint_d, rint_sub)

      loess_preds <- stats::predict(d.loess, nd)
    } else {
      loess_preds <- rep(empty_val, length(nd))
    }

    if(type == "PI") {
      loess_preds <- exp(loess_preds)
    }

    results <- data.frame(predictor = x,
                          date = nd,
                          preds = loess_preds,
                          stringsAsFactors = FALSE)

    # do training data check and replace for gaps, edges, and corners
    results$preds <- apply(results,
                           1,
                           train_check,
                           empty_val = empty_val,
                           train_data = D,
                           e = e)

    return(results)
  }

  if(!all(is.na(st_extent))) {
    if(!is.list(st_extent)) {
      stop("The st_extent argument must be a list object.")
    }
  }

  # load config vars
  e <- load_config(path)

  # static vars and projection info
  PD_MAX_RESOLUTION <- 50
  ll <- "+init=epsg:4326"
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

  # subset centroids
  # need to remove time from st_extent before using to subset
  # but first save the min and max for plotting
  min_date <- as.Date(st_extent$t.min * 366, origin = as.Date('2013-01-01'))
  max_date <- as.Date(st_extent$t.max * 366, origin = as.Date('2013-01-01'))

  st_extent$t.min <- NULL
  st_extent$t.max <- NULL

  tpis_sub <- st_extent_subset(pis, st_extent)
  tpis <- tpis_sub[!is.na(tpis_sub[, 2]), ]
  rm(pis, tpis_sub)

  tpds_sub <- st_extent_subset(pds, st_extent)
  tpds <- tpds_sub[!is.na(tpds_sub$V5), ]
  rm(pds)

  # subset to cover classes
  land.cover.class.codes <- c(1:10, 12, 13, 16)
  lc.tag <- "UMD_FS_C"
  land_covers <- paste(lc.tag, land.cover.class.codes, sep = "")

  water.cover.class.codes <- c(0, 2, 3, 5, 6, 7)
  wc.tag <- "MODISWATER_FS_C"
  water_covers <- paste(wc.tag, water.cover.class.codes, sep = "")

  cover_classes <- c(land_covers, water_covers)
  cover_cols <- e$PREDICTOR_LIST[grep(paste(cover_classes,
                                            collapse="|"),
                                      e$PREDICTOR_LIST)]

  # if only using PLAND and LPI fragtats, subset further
  if(pland_and_lpi_only == TRUE) {
    cover_cols <- cover_cols[grep(pattern = "PLAND|LPI", cover_cols)]
  }

  if(length(cover_cols) == 0) {
    stop("No land or water cover classes to plot.")
  }

  ## Calculate PD slopes for each centroid
  calc_slope <- function(y) {
    x_vals <- as.numeric(tpds_covs[y, (PD_MAX_RESOLUTION+4):(2*PD_MAX_RESOLUTION+3)])
    y_vals <- as.numeric(tpds_covs[y, 3:(PD_MAX_RESOLUTION+2)]) -
      mean(as.numeric(tpds_covs[y, 3:(PD_MAX_RESOLUTION+2)]), na.rm = T)

    if(all(is.na(y_vals))) {
      return(NA)
    }

    nnn <- length(x_vals)
    xxx <- matrix(cbind(rep(1, nnn), y_vals), nnn, 2)
    sl <- .lm.fit(x = xxx, y = x_vals)$coefficients[[2]]

    if(!is.na(sl)) {
      if(sl > 0) {
        return(1)
      } else {
        return(0)
      }
    } else {
      return(NA)
    }
  }

  tpds_covs <- tpds[tpds$V4 %in% cover_cols, ]
  rm(tpds)

  pd_slope <- lapply(X = 1:nrow(tpds_covs), FUN=calc_slope)
  tpds_covs$slope <- unlist(pd_slope)
  rm(pd_slope)

  pd_w_slope <- tpds_covs[, c("stixel.id", "date", "lat", "lon", "stixel_width",
                              "stixel_height", "V4", "slope")]

  # this aggregation is both for speed and stability of loess
  # without this agg, some of the pd trajectories
  # for lesser predictors get a little weird
  pd_mean_slope <- stats::aggregate(pd_w_slope,
                                    by = list(pd_w_slope$V4,
                                              pd_w_slope$date),
                                    FUN = mean,
                                    na.rm = TRUE)
  rm(pd_w_slope)

  # subset and rename
  pd_slopes <- pd_mean_slope[,c("Group.1", "date", "lat", "lon", "stixel_width",
                                "stixel_height", "slope")]
  names(pd_slopes) <- c("predictor", "date", "lat", "lon", "stixel_width",
                        "stixel_height", "slope")
  rm(pd_mean_slope)

  # temporal smoothing of predictor trajectories
  loesses <- lapply(X = unique(pd_slopes$predictor),
                    FUN = loess_fit_and_predict,
                    ext = st_extent,
                    input_data = pd_slopes,
                    type = "PD",
                    e = e)
  smooth_pds <- dplyr::bind_rows(loesses)
  rm(loesses, pd_slopes)

  # categorize positive and negative directionalities
  smooth_pds$direction <- NA
  smooth_pds$direction[smooth_pds$preds >= 0.7] <- 1
  smooth_pds$direction[smooth_pds$preds <= 0.3] <- -1
  smooth_pds$preds <- NULL

  # calculate PIs
  pi_loess <- lapply(X = cover_cols,
                     FUN = loess_fit_and_predict,
                     ext = st_extent,
                     input_data = tpis,
                     type = "PI",
                     e = e)
  smooth_pis <- dplyr::bind_rows(pi_loess)
  rm(pi_loess, tpis)

  # scale PIs
  pi_week_sums <- stats::aggregate(smooth_pis$preds,
                                   by = list(smooth_pis$date),
                                   FUN = sum,
                                   na.rm = TRUE)

  smooth_pis_w_sums <- merge(smooth_pis,
                             pi_week_sums,
                             by.x = "date",
                             by.y = "Group.1")
  rm(smooth_pis, pi_week_sums)

  smooth_pis_w_sums$pi_adj <- smooth_pis_w_sums$preds/smooth_pis_w_sums$x
  smooth_pis_w_sums$x <- NULL
  smooth_pis_w_sums$preds <- NULL

  # mutiply PI x PD
  pipd <- merge(smooth_pis_w_sums, smooth_pds, by=c("predictor", "date"))
  pipd$pidir <- pipd$pi_adj * pipd$direction
  pipd$direction <- NULL
  pipd$pi_adj <- NULL
  rm(smooth_pis_w_sums, smooth_pds)

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

#' Converts cryptic cover class names to readable land cover names
#'
#' Internal function that converts the cryptic predictor class names to
#' readable land cover names.
#'
#' @param cov_names vector; Cover class names to convert.
#' @param by_cover_class logical; Default is FALSE. If TRUE, replaces fragstat
#' cover class name with a name for the cover class as whole.
#' @param pretty logical; Default is FALSE. If TRUE, Converts from capital case
#' to title case.
#'
#' @return A vector of converted names.
#'
#' @keywords internal
#'
#' @examples
#' cns <- c("UMD_FS_C1_1500_PLAND", "MODISWATER_FS_C7_1500_LPI")
#'
#' converted <- convert_classes(cov_names = cns, pretty = TRUE)
convert_classes <- function(cov_names,
                            by_cover_class = FALSE,
                            pretty = FALSE) {

  if(by_cover_class == TRUE) {
    ending <- ""
  } else {
    ending <- "_"
  }

  # subset to cover classes
  land.cover.class.codes <- c(1:10, 12, 13, 16)
  lc.tag <- "UMD_FS_C"
  land_covers <- paste(lc.tag, land.cover.class.codes, ending, sep = "")

  water.cover.class.codes <- c(0, 2, 3, 5, 6, 7)
  wc.tag <- "MODISWATER_FS_C"
  water_covers <- paste(wc.tag, water.cover.class.codes, ending, sep = "")

  both_covers <- c(land_covers, water_covers)

  land_cover_names <- c("EVERGREEN_NEEDLELEAF_FOREST",
                        "EVERGREEN_BROADLEAF_FOREST",
                        "DECIDUOUS_NEEDLELEAF_FOREST",
                        "DECIDUOUS_BROADLEAF_FOREST",
                        "MIXED_FOREST",
                        "CLOSED_SHRUBLANDS",
                        "OPEN_SHRUBLANDS",
                        "WOODY_SAVANNAS",
                        "SAVANNAS",
                        "GRASSLANDS",
                        "CROPLANDS",
                        "URBAN",
                        "BARREN")

  water_cover_names <- c("SHALLOW_OCEAN", "OCEAN_COASTLINES_AND_LAKE_SHORES",
                         "SHALLOW_INLAND_WATER", "DEEP_INLAND_WATER",
                         "MODERATE_OCEAN", "DEEP_OCEAN")

  both_names <- paste(c(land_cover_names, water_cover_names), ending, sep = "")

  converted <- c()
  for(n in cov_names) {
    a <- both_covers[which(!is.na(pmatch(both_covers, n)))]
    b <- both_names[which(!is.na(pmatch(both_covers, n)))]

    conv <- ifelse(length(a) > 0, stringr::str_replace_all(n, a, b), n)
    converted <- c(converted, conv)
  }

  if(pretty == TRUE) {
    converted <- lettercase::str_title_case(
      lettercase::str_lower_case(converted))
  }

  return(converted)
}
