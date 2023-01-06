context("Loading functions")

skip_on_cran()

test_that("get_species_path", {
  p <- get_species_path("example_data")
  expect_true(dir.exists(p))
  expect_error(get_species_path("not-a-real-species"))
})

test_that("load_config", {
  p <- load_config(path)
  expect_is(p, "list")
  expect_true(all(c("bins", "bins_seasonal", "srd_pred_year") %in% names(p)))
  expect_error(load_config("/invalid/path/"))
})


test_that("load_fac_map_parameters", {
  p <- load_fac_map_parameters(path)
  expect_is(p, "list")
  expect_named(p, c("custom_projection", "fa_extent", "res", "fa_extent_sinu",
                    "weekly_bins", "weekly_labels",
                    "seasonal_bins", "seasonal_labels"))

  # check components
  # projection
  expect_is(terra::crs(p$custom_projection), "character")
  # extent
  expect_is(p$fa_extent, "SpatExtent")
  # resolution
  expect_is(p$res, c("numeric", "integer"))
  # sinusoidal extent
  expect_is(p$fa_extent_sinu, "SpatExtent")
  # bins
  expect_is(p$weekly_bins, "numeric")
  expect_is(p$seasonal_bins, "numeric")

  expect_error(load_config("/invalid/path/"))
})
