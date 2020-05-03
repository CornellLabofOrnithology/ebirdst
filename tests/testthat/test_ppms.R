context("PPM functions")

sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.5, 0.6))

context("compute_ppms")

test_that("ebirdst compute_ppms", {
  # expected
  cppms <- compute_ppms(sp_path, ext = lp_extent)

  expect_equal(length(cppms), 3)
  expect_is(cppms, "list")
  expect_named(cppms, c("binary_ppms", "occ_ppms", "abd_ppms"))
  expect_is(cppms$binary_ppms, "data.frame")
  expect_is(cppms$occ_ppms, "data.frame")
  expect_is(cppms$abd_ppms, "data.frame")
  expect_is(cppms$binary_ppms$mc_iteration, "integer")
  expect_is(cppms$occ_ppms$mc_iteration, "integer")
  expect_is(cppms$abd_ppms$mc_iteration, "integer")
  expect_length(cppms$binary_ppms$mc_iteration, 25)
  expect_length(cppms$occ_ppms$mc_iteration, 25)
  expect_length(cppms$abd_ppms$mc_iteration, 25)
})

context("plot_binary_by_time")

# plot binary by time
test_that("ebirdst plot_binary_by_time", {
  # checking metrics
  expect_error(plot_binary_by_time(sp_path, metric = "kappa",
                                   ext = lp_extent,
                                   n_time_periods = 12), NA)

  # wrong metric
  expect_error(plot_binary_by_time(sp_path, metric = "wrongmetric",
                                   ext = lp_extent,
                                   n_time_periods = 12))

  # n_time_periods
  expect_error(plot_binary_by_time(sp_path, metric = "kappa",
                                   ext = lp_extent,
                                   n_time_periods = 1))
})

context("plot_all_ppms")

test_that("ebirdst plot_all_ppms", {
  # expected
  expect_error(plot_all_ppms(path = sp_path, ext = lp_extent), NA)

  # missing temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_error(plot_all_ppms(path = sp_path, ext = nt))

  # broken path
  expect_error(plot_all_ppms(path = "/bad/path/", ext = lp_extent))
})
