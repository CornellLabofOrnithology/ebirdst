#' Load PPM data
#'
#' @export
load_ppm_data <- function(path) {

  # load the test data and assign names
  test_file <- paste(path,
                     "/results/abund_preds/unpeeled_folds/test.pred.ave.txt",
                     sep = "")

  if(!file.exists(test_file)) {
    stop("*_erd.test.data.csv file does not exist in the /data directory.")
  }

  ppm_data <- data.table::fread(test_file)
  ppm_names <- c("data.type", "row.id", "lon", "lat", "date", "obs", "pi.mean",
                 "pi.90", "pi.10", "pi.se", "pi.mu.mean", "pi.mu.90",
                 "pi.mu.10", "pi.mu.se", "pat", "pi.es")
  names(ppm_data) <- ppm_names

  # load ensemble support values that define weekly extent of analysis
  es_dir <- paste(path,
                  "/results/tifs/presentation/",
                  "abundance_ensemble_support_values/",
                  sep = "")
  es_files <- list.files(es_dir)

  eoa_es_data <- data.frame(Week = NA,
                            NorthAmerica = NA,
                            SouthAmerica = NA)

  for (iii in 1:length(es_files)) {
    ttt <- read.csv(paste(es_dir, es_files[iii], sep = ""))
    eoa_es_data[iii, 1] <- as.character(ttt[1, 2])
    eoa_es_data[iii, 2:3] <- ttt[1, 3:4]
  }

  # add day of year
  eoa_es_data$DOY <- as.numeric(format(strptime(x = eoa_es_data$Week,
                                                format = "%m-%d"),
                                       "%j"))

  # smooth the extent of estimate ensemble support values to day

  # init
  pred_DOY <- data.frame(DOY = c(1:366))
  eoa_NorthAmerica <- rep(NA, nrow(pred_DOY))
  eoa_SouthAmerica <- rep(NA, nrow(pred_DOY))

  # check for NA data
  if(sum(!is.na(eoa_es_data$NorthAmerica)) > 10) {
    # Treat missing ES values as the minimum support value
    na_na_replace <- min(eoa_es_data$NorthAmerica, na.rm = TRUE)
    eoa_es_data$NorthAmerica[is.na(eoa_es_data$NorthAmerica)] <- na_na_replace

    # GAM Smooth EOA ES values down to Daily
    s = mgcv::s
    na_gam <- mgcv::gam(NorthAmerica ~ s(DOY, k = 25, bs = "cp", m = 1),
                        gamma = 1.5,
                        data = eoa_es_data,
                        knots = list(DOY = c(1,366)))
    eoa_NorthAmerica <- predict(na_gam, newdata = pred_DOY)
  }

  if(sum(!is.na(eoa_es_data$SouthAmerica)) > 10) {
    # Treat missing ES values as the minimum support value
    sa_na_replace <- min(eoa_es_data$SouthAmerica, na.rm = TRUE)
    eoa_es_data$SouthAmerica[is.na(eoa_es_data$SouthAmerica)] <- sa_na_replace

    # GAM Smooth EOA ES values down to Daily
    s = mgcv::s
    sa_gam <- mgcv::gam(SouthAmerica ~ s(DOY, k = 25, bs = "cp", m = 1),
                        gamma = 1.5,
                        data = eoa_es_data,
                        knots = list(DOY = c(1,366)))
    eoa_SouthAmerica <- predict(sa_gam, newdata = pred_DOY)
  }

  # package and return
  eoa_es_daily <- data.frame(DOY = pred_DOY,
                             NorthAmerica = eoa_NorthAmerica,
                             SouthAmerica = eoa_SouthAmerica  )

  return(list(ppm_data = ppm_data,
              eoa_es_daily =eoa_es_daily))
}
