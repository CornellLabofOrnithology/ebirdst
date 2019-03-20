context("Plotting functions")

sp_path <- download_data("example_data", tifs_only = FALSE, force = TRUE)
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.5, 0.6))

context("plot_pis")
# plot_pis
test_that("ebirdst plot_pis", {
  pis <- load_pis(sp_path)

  # expected with ebirdst_extent
  top_pi <- plot_pis(pis = pis, ext = lp_extent, plot = FALSE)
  expect_is(top_pi, "numeric")

  # check names and lengths work correctly
  expect_true("Mixed Forest ED" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = FALSE,
                     plot = FALSE)
  expect_true("umd_fs_c5_1500_ed" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = FALSE,
                     by_cover_class = TRUE, plot = FALSE)
  expect_true("umd_fs_c5" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = TRUE,
                     by_cover_class = TRUE, plot = FALSE)
  expect_true("Mixed Forest" %in% names(top_pi))

  # checking length
  top_pi <- plot_pis(pis = pis, ext = lp_extent, n_top_pred = 10,
                     plot = FALSE)
  expect_length(top_pi, 10)
  top_pi <- plot_pis(pis = pis, ext = lp_extent, n_top_pred = 10,
                     by_cover_class = TRUE, plot = FALSE)
  expect_length(top_pi, 10)

  # not enough num_top_preds
  expect_error(plot_pis(pis = pis, ext = lp_extent, n_top_pred = 1,
                        plot = FALSE))

  # broken input
  expect_error(plot_pis(pis = pis))
  expect_error(plot_pis(pis = 1:10, ext = lp_extent))
  expect_error(plot_pis(pis = data.frame(), ext = lp_extent))

  # missing temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_error(plot_pis(pis = pis, ext = nt))
})

context("plot_pds")

# plot_pds
test_that("ebirdst plot_pds", {
  pds <- load_pds(sp_path)

  # checking return quality
  pds_data <- plot_pds(pds = pds, predictor = "effort_hrs",
                       ext = lp_extent, plot = FALSE)
  expect_is(pds_data, "data.frame")
  expect_named(pds_data, c("x", "pd_median", "pd_lower", "pd_upper"))
  expect_equal(nrow(pds_data), 50)

  # quantiles working
  pds_data_q <- plot_pds(pds = pds, predictor = "effort_hrs",
                         ext = lp_extent, show_quantiles = TRUE, plot = FALSE)
  expect_is(pds_data_q, "data.frame")
  expect_named(pds_data_q, c("x", "pd_median", "pd_lower", "pd_upper",
                             "pd_lower_quantile", "pd_upper_quantile"))
  expect_equal(nrow(pds_data_q), 50)

  # more params
  expect_is(plot_pds(pds = pds,
                     predictor = "EFFORT_HRS",
                     ext = lp_extent,
                     show_quantiles = TRUE,
                     show_stixel_pds = TRUE), "data.frame")
  expect_is(plot_pds(pds = pds,
                    predictor = "EFFORT_HRS",
                    ext = lp_extent,
                    ylim = c(-1, 1)), "data.frame")
  expect_error(plot_pds(pds = pds,
                      predictor = "EFFORT_HRS",
                      ext = lp_extent,
                      ylim = 1))

  # TODO add more tests for GBM parms

  # missing temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_error(plot_pds(pds = pds, predictor = "effort_hrs",
                        ext = nt, plot = FALSE))
})
