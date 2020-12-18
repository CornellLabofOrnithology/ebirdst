context("Loading functions")

path <- ebirdst_download("example_data", tifs_only = TRUE)

test_that("load_config", {
  l <- load_config(path)
  expect_is(l, "list")
  expect_true(all(c("bins", "bins_seasonal", "SRD_PRED_YEAR") %in% names(l)))
  expect_error(load_config("/invalid/path/"))
})


test_that("load_fac_map_parameters", {
  l <- load_fac_map_parameters(path)
  expect_is(l, "list")
  expect_named(l, c("custom_projection", "fa_extent", "res",
                    "fa_extent_sinu", "abundance_bins"))

  # check components
  # projection
  expect_is(raster::projection(l$custom_projection), "character")
  # extent
  expect_is(l$fa_extent, "Extent")
  # resolution
  expect_is(l$res, "numeric")
  # sinusoidal extent
  expect_is(l$fa_extent_sinu, "Extent")
  # bins
  expect_is(l$abundance_bins, "numeric")

  expect_error(load_config("/invalid/path/"))
})
