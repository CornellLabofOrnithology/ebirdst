context("ebirdst_extent")

test_that("extent works", {
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))

  # spatial
  expect_equal(e$type, "bbox")
  expect_is(e$extent, "bbox")
  expect_equal(as.numeric(e$extent), c(-86, 42, -83, 45))

  # temporal
  expect_is(e$t, "numeric")
  expect_equal(e$t, c(0.5, 0.6))

  # polygon for spatial, dates for temporal
  e <- ebirdst_extent(sf::st_as_sfc(e$extent),
                      t = c("2019-01-01", "2019-12-31"))

  # spatial
  expect_equal(e$type, "polygon")
  expect_is(e$extent, "sfc")

  # temporal
  expect_is(e$t, "numeric")
  expect_equal(round(e$t, 2), c(0, 1))
})

test_that("extent erors", {
  # reversed min max
  expect_error(ebirdst_extent(c(xmin = -83, xmax = -86, ymin = 42, ymax = 45)))
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47, ymax = 45)))

  # missing a corner
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47)))

  # non sequential dates
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47),
                              t = c(1, 0)))

  # invalid dates
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47),
                              t = c(-1, 0)))
})
