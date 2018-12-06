context("Plotting functions")
context("plot_pis")

# plot_pis
test_that("ebirdst plot_pis", {
  pis <- load_pis(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82,
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
  expect_error(plot_pis(path = '~/some/messed/up/path/that/does/not/exist',
                        pis = pis,
                        st_extent = ne_extent),
               "RData file does not exist")

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 45,
                    lon.min = -88,
                    lon.max = -82)
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
               "Minimum latitude is greater than maximum latitude")

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
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "st_extent argument must be a list")
})

context("plot_pds")

# plot_pds
test_that("ebirdst plot_pds", {
  pds <- load_pds(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 45,
                    lon.min = -88,
                    lon.max = -82,
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
               "Minimum latitude is greater than maximum latitude")

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
