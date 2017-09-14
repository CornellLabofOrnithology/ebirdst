

#' Plotting PIs as barplots
#'
#' @export
#' @import ggplot2 mgcv
plot_pis <- function(pis,
                     st_extent,
                     by_cover_class = FALSE,
                     num_top_preds = 50) {

  # subset for exetnt
  ttt <- pis[pis$centroid.date > st_extent$t.min &
               pis$centroid.date <= st_extent$t.max &
               pis$centroid.lat > st_extent$y.min &
               pis$centroid.lat <= st_extent$y.max &
               pis$centroid.lon > st_extent$x.min &
               pis$centroid.lon <= st_extent$x.max, 2:87]

  if(by_cover_class == TRUE) {

    land.cover.class.codes <- c(1:10,12,13,16)
    lc.tag <- "UMD_FS_C"
    water.cover.class.codes <- c(0,2,3,5,6,7)
    wc.tag <- "MODISWATER_FS_C"
    predictor.names <- names(ttt)

    # ---
    ttt.new <- NULL
    for (iii.pred in 1:length(land.cover.class.codes)){
      # iii.pred <- 1
      new.pred.name <- paste(lc.tag,land.cover.class.codes[iii.pred],"_",sep="")
      pred.nindex <- grep(
        new.pred.name,
        x = predictor.names)
      ttt.new <- cbind(ttt.new,
                       apply( ttt[, pred.nindex], 1, mean, na.rm=T))
      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <-
        paste(lc.tag,land.cover.class.codes[iii.pred],sep="")
    }
    for (iii.pred in 1:length(water.cover.class.codes)){
      # iii.pred <- 1
      new.pred.name <- paste(wc.tag,water.cover.class.codes[iii.pred],"_",sep="")
      pred.nindex <- grep(
        new.pred.name,
        x = predictor.names)
      ttt.new <- cbind(ttt.new,
                       apply( ttt[, pred.nindex], 1, mean, na.rm=T))
      ttt.new <- as.data.frame(ttt.new)
      names(ttt.new)[ncol(ttt.new)] <-
        paste(wc.tag,water.cover.class.codes[iii.pred],sep="")
    }
    #head(ttt.new)
    ttt <- ttt.new
  }

  # compute median
  pi_median <- apply(ttt, 2, median, na.rm = T)

  # find the top preds based on function variable num_top_preds
  top_names <- names(pi_median)[order(pi_median,
                                      decreasing = T)][1:num_top_preds]
  top_names <- na.omit(top_names)

  # subset all values based on top_names
  top_pis <- ttt[ ,top_names]

  # munging and filtering for ggplot
  pi_stack <- stack(top_pis)

  # PI's have have spurious large values, NA's and NAN's
  pi_stack$values[!is.numeric(pi_stack$values)] <- NA
  pi_stack <- pi_stack[pi_stack$values < quantile(pi_stack$values,
                                                  probs = c(0.98),
                                                  na.rm = TRUE), ]
  pi_stack <- pi_stack[complete.cases(pi_stack),]

  pi_bars <- ggplot2::ggplot(pi_stack, aes(reorder(ind,
                                                   values,
                                                   FUN=median),
                                           values)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = "Relative PI", x = "")
  pi_bars
}

#' Plot PDs
#'
#' @export
plot_pds <- function(pd_name,
                     pds,
                     st_extent,
                     pointwise_pi = FALSE,
                     stixel_pds = TRUE,
                     k.cont.res = 25,
                     gbm.n.trees = 500,
                     nnn.bs = 100,
                     equivalent.ensemble.ss = 10,
                     ci.alpha = 0.05,
                     mean.all.data = FALSE) {

  PD_MAX_RESOLUTION <- 50

  # ----------------------
  # PD_NAME <- "EFFORT_HRS"
  # PD_NAME <- "UMD_FS_C4_1500_PD"
  # pipd.data.list,
  # st.extent.list = st.extent.list.E
  # pointwise.ci = T
  # k.cont.res <- 25
  # -------------------
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
  # ----------------------

  pd.index <- pds$centroid.date > st_extent$t.min &
    pds$centroid.date <= st_extent$t.max &
    pds$centroid.lat > st_extent$y.min &
    pds$centroid.lat <= st_extent$y.max &
    pds$centroid.lon > st_extent$x.min &
    pds$centroid.lon <= st_extent$x.max

  pd_vec <- pds[ pd.index, 	]
  # Select PD Variable
  var_pd <- pd_vec[pd_vec$V4 == pd_name,]
  # Clean
  var_pd <- var_pd[!is.na(var_pd$V5), ]
  # Each Column is one replicate estimate of PD
  # 	x = x coordinate values
  # 	y = y coordinate values
  pd.x <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.y <- matrix(NA, PD_MAX_RESOLUTION, nrow(var_pd))
  pd.mean <- rep(NA, nrow(var_pd))
  for (rid in 1:nrow(var_pd)) {
    #rid <- 100
    pd.x[,rid] <- as.numeric(
      var_pd[rid, (PD_MAX_RESOLUTION+4):(2*PD_MAX_RESOLUTION+3)] )
    ttt <- as.numeric(var_pd[rid, 3:(PD_MAX_RESOLUTION+2)])
    pd.mean[rid] <- mean(ttt, na.rm=T)
    pd.y[,rid] <- ttt - pd.mean[rid]
  }
  pd.x <- as.data.frame(pd.x)
  pd.y <- as.data.frame(pd.y)
  # Compute Prediction Design for 1D PD
  ttt <- data.frame(
    x = stack(pd.x)[,1],
    y = stack(pd.y)[,1] )
  nd <- data.frame( x = seq(
    from = quantile(
      ttt$x, probs = x.tail.level, na.rm=T),
    to = quantile(
      ttt$x, probs = 1 - x.tail.level, na.rm=T),
    length = nd.pred.size ) )
  # PLOT STIXEL PD Replicates or just set up plot
  if (stixel_pds){
    matplot(
      jitter(as.matrix(pd.x), amount=0.00),
      pd.y,
      xlab = pd_name,
      ylab = "Deviation E(Logit Occurrence)",
      type= "l",
      #pch = 15,
      #cex = 2.0,
      lwd = 5,
      lty = 1,
      col=alpha("black", .025))
    #ylim = quantile(pd.y, probs=c(0.01, 0.99), na.rm=T))
  }
  if (!stixel_pds){
    plot(
      pd.x[,1],
      pd.y[,1],
      xlab = pd_name,
      ylab = "Deviation E(Logit Occurrence)",
      type = "n")
    #ylim = quantile(pd.y, probs=c(0.01, 0.99), na.rm=T))
  }
  abline(0,0, col="black", lwd=2)
  # -----------------
  # GBM Qunatiles
  # -----------------
  if(pointwise_pi) {
    ttt <- data.frame(
      x = stack(pd.x)[,1],
      y = stack(pd.y)[,1] )
    d.ul <- gbm::gbm(
      y ~ x,
      data = ttt,
      distribution =
        list(name="quantile",alpha=(1-gbm.tail.prob)),
      n.trees = gbm.n.trees,
      interaction.depth = 4,
      shrinkage = 0.05,
      bag.fraction = 0.5,
      train.fraction = 1.0,
      cv.folds = gbm.cv.folds,
      verbose=F,
      n.cores = 1)
    #best.iter <- gbm.perf(d.ul, method="cv", plot.it=F)
    #print(best.iter)
    d.ll <- gbm::gbm(
      y ~ x,
      data = ttt,
      distribution =
        list(name="quantile",alpha=gbm.tail.prob),
      n.trees = gbm.n.trees,
      interaction.depth = 4,
      shrinkage = 0.05,
      bag.fraction = 0.5,
      train.fraction = 1.0,
      cv.folds = gbm.cv.folds,
      verbose=F,
      n.cores = 1)
    #best.iter <- gbm.perf(d.ll, method="cv", plot.it=F)
    #cat("LL:", jjj, best.iter,"\n")
    # ------------
    #best.iter <- 200
    t.ul <- predict(d.ul,
                    newdata= nd,
                    n.trees=best.iter)
    t.ll <- predict(d.ll,
                    newdata= nd,
                    n.trees=best.iter)
    # ------------
    poly.x <- c(nd[,1], rev(nd[,1]))
    poly.y <- c(t.ll, rev(t.ul))
    polygon( poly.x, poly.y, col=alpha("red", .25), border=F)
  } # END if (pointwise.ci){

  # -----------------
  # GAM Pointwise CI for conditional mean estimate
  # via bootstrapping
  # -----------------
  if (pointwise_pi){
    # nnn.bs <- 25
    # nd.pred.size <- 50
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
      random.index <-  matrix(
        (rbinom(
          n = nrow(pd.x)*ncol(pd.x),
          size = 1,
          prob = equivalent.ensemble.ss*PD_MAX_RESOLUTION/
            nrow(pd.x)/ncol(pd.x)) == 1),
        nrow(pd.x),
        ncol(pd.x))
      # dim(random.index)
      # head(random.index)
      # sum(random.index)
      ttt <- data.frame(
        x = pd.x[random.index],
        y = pd.y[random.index] )
      #x = stack(pd.x[row.index, bs.index])[,1],
      #y = stack(pd.y[row.index, bs.index])[,1] )

      s = mgcv::s
      d.gam <- mgcv::gam(
        y ~ s(x, k = k.cont.res, bs="ds", m=1),
        data = ttt,
        gamma = 1.5 )
      # summary(d.gam)
      bs.gam.pred[, iii.bs] <- predict(d.gam, newdata = nd, se=F)
      #lines(nd[,1],bs.gam.pred[, iii.bs], col=alpha("green", 0.5), lwd=2 )
    }
    t.ul <- apply(bs.gam.pred, 1, quantile, probs= 1-ci.alpha, na.rm=T)
    t.ll <- apply(bs.gam.pred, 1, quantile, probs= ci.alpha, na.rm=T)
    t.median <- apply(bs.gam.pred, 1, quantile, probs= 0.5, na.rm=T)
    # ------------
    poly.x <- c(nd[,1], rev(nd[,1]))
    poly.y <- c(t.ll, rev(t.ul))
    polygon( poly.x, poly.y, col=alpha("blue", .25), border=F)
    lines(nd[,1], t.median,
          col=alpha("darkorange", 1.0), lwd=2 )
  }
  # GAM CONDITIONAL MEAN - ALL DATA
  if (mean.all.data){
    ttt <- data.frame(
      x = stack(pd.x)[,1],
      y = stack(pd.y)[,1] )
    d.gam <- mgcv::gam(
      y ~ s(x, k = k.cont.res, bs="ds", m=1),
      data = ttt,
      gamma = 1.5 )
    # summary(d.gam)
    p.gam <- predict(d.gam,
                     newdata = nd,
                     se=T)
    polygon(
      x = c(nd[,1], rev(nd[,1])),
      y = c( p.gam$fit + 2*p.gam$se.fit,
             rev(p.gam$fit - 2*p.gam$se.fit) ),
      col = alpha("lightblue", 0.25),
      border = NA)
    lines( nd[,1], p.gam$fit, col = alpha("yellow", 0.75), lwd=3)
  }
  return( list(
    t.median = t.median,
    t.ul = t.ul,
    t.ll = t.ll))
}

# ' Function used by cake_plot()
# '
train_check <- function(x) {
  train_index <- train_data.date$train_data.date >= (as.numeric(x["date"])-full_week) &
    train_data.date$train_data.date <= (as.numeric(x["date"])+full_week)

  if(all(train_index == FALSE)) {
    return(empty_val)
  } else {
    min_train <- min(train_data.date$train_data.date[train_index], na.rm=TRUE)
    max_train <- max(train_data.date$train_data.date[train_index], na.rm=TRUE)

    # end of year present problems for checking to see if there's data
    # on both sides of the prediction date, this adjusts for that
    pred_date <- as.numeric(x["date"])

    if(pred_date == 0.01) {
      pred_date <- 0.02
    } else if(pred_date == 0.99) {
      pred_date <- 0.98
    }

    if( min_train < pred_date & max_train > pred_date )  {
      # if there is data on both sides, return prediction
      return(as.numeric(x["preds"]))
    } else {
      # otherwise return log "zeros"
      return(empty_val)
    }

  }
}

# ' Make a cake plot of the pis and pds
# '
# ' @export
cake_plot <- function(path,
                      pis,
                      pds,
                      st_extent,
                      by_cover_class = FALSE) {

  # load config vars
  e <- load_config(path)

  # subset centroids
  tpis <- pis[pis$centroid.lat > st_extent$y.min &
                pis$centroid.lat <= st_extent$y.max &
                pis$centroid.lon > st_extent$x.min &
                pis$centroid.lon <= st_extent$x.max &
                !is.na(pis[,2]), ]

  tpds <- pds[pds$centroid.lat > st_extent$y.min &
                pds$centroid.lat <= st_extent$y.max &
                pds$centroid.lon > st_extent$x.min &
                pds$centroid.lon <= st_extent$x.max &
                !is.na(pds$V5), ]

  # subset to cover classes
  land.cover.class.codes <- c(1:10,12,13,16)
  lc.tag <- "UMD_FS_C"
  land_covers <- paste(lc.tag, land.cover.class.codes, sep = "")

  water.cover.class.codes <- c(0,2,3,5,6,7)
  wc.tag <- "MODISWATER_FS_C"
  water_covers <- paste(wc.tag, water.cover.class.codes, sep = "")

  cover_classes <- c(land_covers, water_covers)

  cover_cols <- e$PREDICTOR_LIST[grep(paste(cover_classes,
                                            collapse="|"),
                                      e$PREDICTOR_LIST)]

  # need to extract directionality
  # area weighted loess?
  # minimum number of points?
  PD_MAX_RESOLUTION <- 50

  ## Calculate slopes for each
  #### PD
  calc_slope <- function(y) {
    x <- as.data.frame(tpds[y, ])

    x_vals <- as.numeric(x[(PD_MAX_RESOLUTION+4):(2*PD_MAX_RESOLUTION+3)])
    y_vals <- as.numeric(x[3:(PD_MAX_RESOLUTION+2)]) -
      mean(as.numeric(x[3:(PD_MAX_RESOLUTION+2)]), na.rm = T)

    sm <- lm(x_vals ~ y_vals)

    sl <- sm$coefficients[[2]]

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

  pd_slope <- lapply(X = 1:nrow(tpds), FUN=calc_slope)
  tpds$slope <- unlist(pd_slope)
  rm(pd_slope)

  pd_w_slope <- tpds[tpds$V4 %in% cover_cols, c("stixel.id",
                                                "centroid.date",
                                                "V4",
                                                "slope")]
  pd_mean_slope <- aggregate(pd_w_slope,
                             by = list(pd_w_slope$V4, pd_w_slope$centroid.date),
                             FUN = mean,
                             na.rm = TRUE)

  pd_slopes <- pd_mean_slope[,c("Group.1", "centroid.date", "slope")]
  names(pd_slopes) <- c("predictor", "date", "slope")

  # need to loess fit and predict for each variable for uniform date set
  SRD_DATE_VEC <- seq(from = 0, to= 1, length= 52 +1)
  SRD_DATE_VEC <- (SRD_DATE_VEC[1:52] + SRD_DATE_VEC[2:(52+1)])/2
  SRD_DATE_VEC <- round(SRD_DATE_VEC,digits=2)

  nd <- data.frame(date=SRD_DATE_VEC)

  fit_and_predict_pd <- function(x) {
    D <- pd_slopes[pd_slopes$predictor == x, ]
    D$predictor <- NULL

    d.loess <- loess(formula = "slope ~ date",
                     defree = 1,
                     data = D)

    loess_preds <- predict(d.loess, nd)


    # do training data check



    results <- data.frame(predictor = x,
                          date = nd,
                          smooth_slopes = loess_preds)



    return(results)
  }

  loesses <- lapply(X=unique(pd_slopes$predictor), FUN=fit_and_predict_pd)

  smooth_pds <- dplyr::bind_rows(loesses)

  smooth_pds$direction <- NA
  smooth_pds$direction[smooth_pds$smooth_slopes >= 0.7] <- 1
  smooth_pds$direction[smooth_pds$smooth_slopes <= 0.3] <- -1
  smooth_pds$smooth_slopes <- NULL

  # aggregate by stixel.id and $V4 calculating mean with na.rm

  # calculate mean PI
  pi_means <- as.data.frame(colMeans(tpis[, cover_cols], na.rm = TRUE))

  fit_and_predict_pi <- function(x) {
    D <- tpis[,c(x, "centroid.date")]
    names(D) <- c("predictor", "date")
    D$predictor <- log(D$predictor + 0.001)

    d.loess <- loess(formula = "predictor ~ date",
                     defree = 1,
                     data = D)

    loess_preds <- predict(d.loess, nd)

    results <- data.frame(predictor = x,
                          date = nd,
                          smooth_pis = exp(loess_preds),
                          stringsAsFactors = FALSE)

    return(results)
  }

  pi_loess <- lapply(X=cover_cols, FUN=fit_and_predict_pi)
  smooth_pis <- dplyr::bind_rows(pi_loess)

  # scale PIs
  pi_week_sums <- aggregate(smooth_pis$smooth_pis,
                            by = list(smooth_pis$date),
                            FUN = sum,
                            na.rm = TRUE)

  smooth_pis_w_sums <- merge(smooth_pis,
                             pi_week_sums,
                             by.x = "date",
                             by.y="Group.1")

  smooth_pis_w_sums$pi_adj <- smooth_pis_w_sums$smooth_pis/smooth_pis_w_sums$x
  smooth_pis_w_sums$x <- NULL
  smooth_pis_w_sums$smooth_pis <- NULL


  # mutiply PI x PD
  pipd <- merge(smooth_pis_w_sums, smooth_pds, by=c("predictor", "date"))
  pipd$pidir <- pipd$pi_adj * pipd$direction
  pipd$direction <- NULL
  pipd$pi_adj <- NULL

  if(by_cover_class == TRUE) {
    # add column of class

    pipd <- pipd[grep("*PLAND|*LPI", pipd$predictor),]


    return_class <- function(x) {
      y <- x["predictor"]

      spls <- strsplit(y, "_")

      return(paste(spls[[1]][1], spls[[1]][2], spls[[1]][3], sep="_"))
    }


    pipd$class <- apply(pipd, 1, return_class)

    # sum pidir by class

    pipd_agg <- aggregate(pipd$pidir,
                          by = list(pipd$class, pipd$date),
                          FUN = sum,
                          na.rm = TRUE)
    names(pipd_agg) <- c("predictor", "date", "pidir")

    pipd <- pipd_agg

  }

  # calculate absolute maxes of
  absmax <- aggregate(pipd$pidir,
                      by = list(pipd$predictor),
                      FUN = function(x) { max(abs(x), na.rm = TRUE)})

  short_set <- absmax[!is.infinite(absmax$x) & absmax$x > 0.01,]
  pipd_short <- pipd[pipd$predictor %in% short_set$Group.1, ]

  pipd_short <- pipd

  # sum by cover class (or not?)

  # agg colors
  agg_colors <- c("blue", "chartreuse", "cyan", "cyan3", "blue2", "blue4",
                  "#48D8AE", "#FCF050", "#E28373", "#C289F3", "#E6E6E6",
                  "#19AB81", "#88DA4A", "#5AB01A", "#A3B36B", "#D1BB7B",
                  "#E1D4AC", "#DCA453", "#EACA57")

  unagg_colors <- c("blue", "blue", "blue", "blue", "chartreuse", "chartreuse", "chartreuse", "chartreuse", "cyan", "cyan", "cyan", "cyan3",  "cyan3",  "cyan3",  "cyan3", "blue2", "blue2", "blue2", "blue2", "blue4", "blue4", "blue4", "blue4",
                    "#48D8AE", "#48D8AE", "#FCF050", "#FCF050", "#FCF050", "#FCF050", "#E28373", "#E28373", "#E28373", "#E28373", "#C289F3", "#C289F3",  "#C289F3", "#C289F3", "#E6E6E6", "#E6E6E6", "#E6E6E6", "#E6E6E6",
                    "#19AB81", "#19AB81", "#19AB81", "#19AB81", "#88DA4A", "#88DA4A", "#88DA4A", "#88DA4A", "#5AB01A", "#5AB01A", "#5AB01A", "#5AB01A", "#A3B36B", "#A3B36B", "#A3B36B", "#A3B36B", "#D1BB7B", "#D1BB7B",
                    "#E1D4AC", "#E1D4AC", "#E1D4AC", "#E1D4AC", "#DCA453", "#DCA453", "#DCA453", "#DCA453", "#EACA57", "#EACA57", "#EACA57", "#EACA57" )

  # ggplot
  wave <- ggplot2::ggplot(pipd_short, ggplot2::aes(x=date,
                                                   y=pidir,
                                                   group=predictor,
                                                   fill=predictor)) +
    ggplot2::geom_area() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::scale_fill_manual(values=agg_colors)
    #ggplot2::theme(legend.position = "none")

  wave
}
