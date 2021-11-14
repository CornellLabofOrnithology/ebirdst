context("Stixel functions")

skip_on_cran()

path <- ebirdst_download("example_data", tifs_only = FALSE)
e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                    t = c(0.5, 0.6))
stx <- load_stixels(path)
stx_sf <- load_stixels(path, return_sf = TRUE)

test_that("load_stixels", {
  # expectations
  expect_is(stx, "data.frame")
  expect_equal(ncol(stx), 96)
  expect_gt(nrow(stx), 0)
  expect_true(all(c("stixel_id", "lon", "lat", "date") %in% names(stx)))

  # spatial format
  expect_is(stx_sf, "sf")
  expect_is(sf::st_geometry(stx_sf), "sfc_POINT")
  expect_equal(ncol(stx_sf), 95)
  expect_gt(nrow(stx_sf), 0)
  expect_true(all(c("stixel_id", "date") %in% names(stx_sf)))

  # invalid input
  expect_error(load_stixels("/invalid/path"))
})

test_that("stixelize", {
  p <- stixelize(stx)
  p_sf <- stixelize(stx_sf)

  # expectations
  expect_is(p, "sf")
  expect_is(sf::st_geometry(p), "sfc_POLYGON")
  expect_equal(ncol(stx) + 1, ncol(p))
  expect_equal(nrow(stx), nrow(p))
  expect_true(all(c("stixel_id", "lon", "lat", "date") %in% names(p)))

  # spatial format
  expect_is(p_sf, "sf")
  expect_is(sf::st_geometry(p_sf), "sfc_POLYGON")
  expect_equal(ncol(stx) + 1, ncol(p_sf))
  expect_equal(nrow(stx), nrow(p_sf))
  expect_true(all(c("stixel_id", "lon", "lat", "date") %in% names(p_sf)))

  # invalid input
  expect_error(stixelize(dplyr::select(stx, -lat)))
})

test_that("stixel_footprint", {
  fp <- stixel_footprint(path, ext = e)

  # expectations
  expect_is(fp, "stixel_footprint")
  expect_is(fp$footprint, "RasterLayer")
  expect_is(fp$centroids, "sf")
  expect_is(sf::st_geometry(fp$centroids), "sfc_POINT")
  expect_is(fp$extent, "ebirdst_extent")

  # invalid input
  expect_error(stixel_footprint("/invalid/path"))
})

test_that("plot stixel_footprint", {
  fp <- stixel_footprint(path, ext = e)
  expect_silent(plot(fp))
})
