context("PPM functions")

test_that("stemhelper compute_ppms", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  expect_error(compute_ppms(sp_path), NA)

  cppms <- compute_ppms(sp_path)

  expect_equal(length(cppms), 3)
  expect_is(cppms, "list")
  expect_is(cppms$binary_stats, "data.frame")
  expect_is(cppms$occ_stats, "data.frame")
  expect_is(cppms$count_stats, "data.frame")
  expect_is(cppms$binary_stats$mc, "integer")
  expect_is(cppms$occ_stats$mc, "integer")
  expect_is(cppms$count_stats$mc, "integer")
  expect_equal(length(cppms$binary_stats$mc), 25)
  expect_equal(length(cppms$occ_stats$mc), 25)
  expect_equal(length(cppms$count_stats$mc), 25)

  # with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent), NA)

  cppms <- compute_ppms(sp_path, st_extent = ne_extent)

  expect_equal(length(cppms), 3)
  expect_is(cppms, "list")
  expect_is(cppms$binary_stats, "data.frame")
  expect_is(cppms$occ_stats, "data.frame")
  expect_is(cppms$count_stats, "data.frame")
  expect_is(cppms$binary_stats$mc, "integer")
  expect_is(cppms$occ_stats$mc, "integer")
  expect_is(cppms$count_stats$mc, "integer")
  expect_equal(length(cppms$binary_stats$mc), 25)
  expect_equal(length(cppms$occ_stats$mc), 25)
  expect_equal(length(cppms$count_stats$mc), 25)

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent), NA)

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "st_extent argument must be a list")
})

# plot binary by time
test_that("stemhelper plot_binary_by_time", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  # checking metrics
  expect_error(plot_binary_by_time(sp_path, metric = "Kappa"), NA)
  expect_error(plot_binary_by_time(sp_path, metric = "AUC"), NA)
  expect_error(plot_binary_by_time(sp_path, metric = "Sensitivity"), NA)
  expect_error(plot_binary_by_time(sp_path, metric = "Specificity"), NA)

  # wrong metric
  expect_error(plot_binary_by_time(sp_path, metric = "WrongMetric"),
               "Predictive performance metric must be one of")

  # with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)

  expect_error(plot_binary_by_time(path = sp_path,
                                   metric = "Kappa",
                                   st_extent = ne_extent), NA)

  # n_time_periods
  expect_error(plot_binary_by_time(path = sp_path,
                                   metric = "Kappa",
                                   st_extent = ne_extent,
                                   n_time_periods = 1),
               "n_time_periods argument must be more than 1")
})

# plot all ppms

test_that("stemhelper plot_all_ppms", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)

  # expected
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent), NA)

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent),
               "Must provide temporal limits")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent),
               "st_extent argument must be a list")

  # broken path

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(plot_all_ppms(path = sp_path, st_extent = ne_extent),
               "file does not exist")
})
