context("Loading and subsetting raster data")

skip_on_cran()

path <- ebirdst_download("example_data", tifs_only = TRUE)

test_that("load_raster", {
  # weekly
  abd <- load_raster(path, "abundance")
  expect_is(abd, "RasterStack")
  expect_equal(raster::nlayers(abd), 52)

  # seasonal
  abd <- load_raster(path, "abundance_seasonal")
  expect_is(abd, "RasterStack")
  expect_equal(raster::nlayers(abd), 4)
  expect_named(abd, c("breeding", "postbreeding_migration",
                      "nonbreeding", "prebreeding_migration"))

  # template
  t <- load_raster(path, "template")
  expect_is(t, "RasterLayer")
  expect_equal(raster::nlayers(t), 1)

  # resolutions
  hr <- load_raster(path, "abundance", resolution = "hr")
  mr <- load_raster(path, "abundance", resolution = "mr")
  lr <- load_raster(path, "abundance", resolution = "lr")
  expect_equal(round(raster::res(hr) * 3), round(raster::res(mr)))
  expect_equal(round(raster::res(hr) * 9), round(raster::res(lr)))
})

test_that("load_raster error", {
  # weekly
  expect_error(load_raster(path, "error"))
  expect_error(load_raster("/invalid/path/", "aundance"))
})

test_that("label_raster_stack", {
  abd <- load_raster(path, "abundance")

  # expected
  expect_match(names(abd), "^w[0-9]{4}\\.[0-9]{2}\\.[0-9]{2}")

  # error
  expect_error(label_raster_stack(abd[[1:3]]))
  expect_error(label_raster_stack(abd[[1]]))

  # expected
  date_vector <- parse_raster_dates(abd)
  expect_is(date_vector, "Date")
})

test_that("subset raster", {
  abd <- load_raster(path, "abundance")
  e <- ebirdst_extent(c(xmin = -86, xmax = -83,ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))
  abd_sub <- ebirdst_subset(abd, e)

  # expected
  expect_gt(raster::ncell(abd_sub), 1)
  expect_equal(raster::nlayers(abd_sub), 5)
  expect_is(abd_sub, "RasterBrick")

  # without temporal info
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  abd_sub <- ebirdst_subset(abd, e)
  # expected
  expect_gt(raster::ncell(abd_sub), 1)
  expect_equal(raster::nlayers(abd_sub), 52)
  expect_is(abd_sub, "RasterBrick")

  # extent with polygon
  e_poly <- ebirdst_extent(sf::st_as_sfc(e$extent), t = c(0.5, 0.6))
  suppressWarnings(abd_sub <- ebirdst_subset(abd, e_poly))
  expect_gt(raster::ncell(abd_sub), 1)
  expect_equal(raster::nlayers(abd_sub), 5)
  expect_is(abd_sub, "RasterBrick")
})
