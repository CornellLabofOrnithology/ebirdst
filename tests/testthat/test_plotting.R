context("Plotting functions")

# plot_pis
test_that("stemhelper plot_pis", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pis <- load_pis(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        by_cover_class = TRUE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 10),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = TRUE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        pretty_names = FALSE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = TRUE,
                        print_plot = FALSE),
               NA)

  # checking return quality
  expect_equal(length(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 25,
                        return_top = TRUE)),
               25)
  expect_is(plot_pis(path = sp_path,
                     pis = pis,
                     st_extent = ne_extent,
                     num_top_preds = 25,
                     return_top = TRUE),
               "character")

  # nothing to do
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = FALSE,
                        print_plot = FALSE),
               "Both print and return params are FALSE")

  # not enough num_top_preds
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 1),
               "num_top_preds must be greater than 1")

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "RData file does not exist")

  sp_path <- paste(root_path, species, sep = "")

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(plot_pis(path = sp_path,
                     pis = pis,
                     st_extent = ne_extent),
               "Must provide t.min and t.max")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(stack_stem(sp_path,
                          variable = "abundance_umean",
                          st_extent = ne_extent),
               "st_extent argument must be a list")
})

# plot_pds
test_that("stemhelper plot_pds", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pds <- load_pds(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent), NA)

  # checking return quality
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent), "list")
  expect_equal(length(plot_pds(pd_name = "EFFORT_HRS",
                               pds = pds,
                               st_extent = ne_extent)$pointwise), 3)
  expect_equal(length(plot_pds(pd_name = "EFFORT_HRS",
                               pds = pds,
                               st_extent = ne_extent)$pointwise$t.median), 50)
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent)$pointwise, "list")
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent)$pointwise$t.median, "numeric")

  # params
  # quantiles
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        plot_quantiles = TRUE), NA)
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent,
                     plot_quantiles = TRUE), "list")
  expect_equal(length(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent,
                     plot_quantiles = TRUE)), 2)
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent,
                     plot_quantiles = TRUE)$quantiles, "list")
  expect_equal(length(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent,
                     plot_quantiles = TRUE)$quantiles), 2)
  expect_is(plot_pds(pd_name = "EFFORT_HRS",
                     pds = pds,
                     st_extent = ne_extent,
                     plot_quantiles = TRUE)$quantiles$t.ul, "numeric")
  expect_equal(length(plot_pds(pd_name = "EFFORT_HRS",
                               pds = pds,
                               st_extent = ne_extent,
                               plot_quantiles = TRUE)$quantiles$t.ul), 50)

  # more params
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        plot_quantiles = TRUE,
                        pointwise_pi = FALSE), NA)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        stixel_pds = TRUE), NA)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        ylim = 1), "invalid 'ylim' value")
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        print_plot = FALSE), NA)

  # TODO add more tests for GBM parms

  # nothing to do
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent,
                        pointwise_pi = FALSE),
               "Nothing to plot")

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent),
               "Must provide t.min and t.max")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(plot_pds(pd_name = "EFFORT_HRS",
                        pds = pds,
                        st_extent = ne_extent),
               "st_extent argument must be a list")
})

# cake_plots
test_that("stemhelper cake_plot", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pis <- load_pis(sp_path)
  pds <- load_pds(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent), NA)

  # params
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent,
                         by_cover_class = FALSE), NA)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent,
                         pland_and_lpi_only = FALSE), NA)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent,
                         by_cover_class = FALSE,
                         pland_and_lpi_only = FALSE), NA)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent,
                         return_data = TRUE), NA)

  # checking return data
  cp <- cake_plot(path = sp_path,
                  pis = pis,
                  pds = pds,
                  st_extent = ne_extent,
                  return_data = TRUE)

  expect_is(cp, "data.frame")
  expect_gt(nrow(cp), 0)
  expect_is(cp$predictor, "character")
  expect_is(cp$date, "numeric")
  expect_is(cp$pidir, "numeric")
  expect_is(cp$labels, "character")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(cake_plot(path = sp_path,
                         pis = pis,
                         pds = pds,
                         st_extent = ne_extent),
               "st_extent argument must be a list")

})
