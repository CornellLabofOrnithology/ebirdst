context("Loading and subsetting raster data")

skip_on_cran()

test_that("load_raster", {
  # weekly
  abd <- load_raster(path, "abundance", resolution = "lr")
  expect_is(abd, "SpatRaster")
  expect_equal(terra::nlyr(abd), 52)

  # check labellling
  expect_match(names(abd), "^[0-9]{4}-[0-9]{2}-[0-9]{2}")
  expect_is(parse_raster_dates(abd))

  # seasonal
  abd <- load_raster(path, "abundance", period = "seasonal", resolution = "lr")
  expect_is(abd, "SpatRaster")
  expect_equal(terra::nlyr(abd), 4)
  expect_named(abd, c("breeding", "nonbreeding",
                      "prebreeding_migration", "postbreeding_migration"))
})

test_that("load_raster error", {
  # weekly
  expect_error(load_raster(path, "error"))
  expect_error(load_raster("/invalid/path/", "aundance"))
})

test_that("subset raster", {
  abd <- load_raster(path, "abundance", resolution = "lr")
  e <- ebirdst_extent(c(xmin = -86, xmax = -83,ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))
  abd_sub <- ebirdst_subset(abd, e)

  # expected
  expect_gt(terra::ncell(abd_sub), 1)
  expect_equal(terra::nlyr(abd_sub), 5)
  expect_is(abd_sub, "SpatRaster")

  # without temporal info
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  abd_sub <- ebirdst_subset(abd, e)
  # expected
  expect_gt(terra::ncell(abd_sub), 1)
  expect_equal(terra::nlyr(abd_sub), 52)
  expect_is(abd_sub, "SpatRaster")

  # extent with polygon
  e_poly <- ebirdst_extent(sf::st_as_sfc(e$extent), t = c(0.5, 0.6))
  suppressWarnings(abd_sub <- ebirdst_subset(abd, e_poly))
  expect_gt(terra::ncell(abd_sub), 1)
  expect_equal(terra::nlyr(abd_sub), 5)
  expect_is(abd_sub, "SpatRaster")
})
