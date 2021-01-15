context("Plotting functions")
skip_on_cran()
options(warn = -1)

sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
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
  expect_true("Deciduous Broadleaf Forests PLAND" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = FALSE,
                     plot = FALSE)
  expect_true("mcd12q1_lccs1_fs_c14_1500_pland" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = FALSE,
                     by_cover_class = TRUE, plot = FALSE)
  expect_true("mcd12q1_lccs1_fs_c14" %in% names(top_pi))
  top_pi <- plot_pis(pis = pis, ext = lp_extent, pretty_names = TRUE,
                     by_cover_class = TRUE, plot = FALSE)
  expect_true("Deciduous Broadleaf Forests" %in% names(top_pi))

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
