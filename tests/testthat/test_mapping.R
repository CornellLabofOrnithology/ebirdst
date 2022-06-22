context("Mapping functions")

skip_on_cran()

e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                    t = c(0.5, 0.6))

test_that("calc_full_extent", {
  abd <- load_raster(path, "abundance", resolution = "lr")
  abd <- ebirdst_subset(abd, e)
  fe <- calc_full_extent(abd)

  # expectations
  expect_is(calc_full_extent(abd), "Extent")
  expect_is(calc_full_extent(abd[[2]]), "Extent")

  # not a raster object
  bb <- sf::st_bbox(abd)
  expect_error(calc_full_extent(bb))
})

test_that("calc_bins", {
  abd <- load_raster(path, "abundance", resolution = "lr")
  cnt <- load_raster(path, "count", resolution = "lr")
  abd <- ebirdst_subset(abd, e)
  cnt <- ebirdst_subset(cnt, e)

  # expectations
  bins <- calc_bins(abd, cnt)
  expect_is(bins, "numeric")
  expect(all(bins > 0), "bins must be positive")

  # single raster layer
  bins <- calc_bins(abd[[1]], cnt[[1]])
  expect_is(bins, "numeric")
  expect(all(bins > 0), "bins must be positive")

  # invalid input
  abd_na <- abd[[2]]
  suppressWarnings(abd_na[] <- NA)
  expect_error(calc_bins(abd_na, cnt[[2]]))
  expect_error(calc_bins(abd_na))
})

test_that("abundance_palette", {
  expect_is(abundance_palette(n = 10), "character")
  expect_length(abundance_palette(n = 10), 10)
  expect_length(abundance_palette(n = 10, season = "breeding"), 10)
  expect_match(abundance_palette(n = 10, season = "nonbreeding"),
               "#[0-9A-F]{6}")
})
