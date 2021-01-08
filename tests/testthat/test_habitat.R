context("Habitat associations")

skip_on_cran()
skip_on_appveyor()
skip_on_travis()

path <- ebirdst_download("example_data", tifs_only = FALSE)
e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43))
habitat <- ebirdst_habitat(path, ext = e, n_predictors = 10)

test_that("ebirdst_habitat", {
  # expected
  expect_is(habitat, "ebirdst_habitat")
  expect_is(habitat, "data.frame")
  expect_named(habitat, c("predictor", "date", "importance", "direction"))
  expect_is(habitat$predictor, "character")
  expect_is(habitat$date, "numeric")
  expect_is(habitat$importance, "numeric")
  expect_is(habitat$direction, "numeric")
  expect_length(unique(habitat$predictor), 10)
  expect_length(unique(habitat$date), 52)

  # invalid inputs
  expect_error(ebirdst_habitat("/invalid/path", ext = e))
})

test_that("ebirdst_habitat extent", {
  # temporal extent ignored
  e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43),
                      t = c(0.5, 0.75))
  h2 <- ebirdst_habitat(path, ext = e, n_predictors = 10)
  expect_identical(habitat, h2)

  # spatial extent required
  expect_error(ebirdst_habitat(path))
})

test_that("ebirdst_habitat n_predictors", {
  h15 <- ebirdst_habitat(path, ext = e, n_predictors = 15)
  expect_length(unique(h15$predictor), 15)

  suppressWarnings({
    h100 <- ebirdst_habitat(path, ext = e, n_predictors = 100)
    expect_lte(length(unique(h100$predictor)), 100)
  })

  # n_predictors must be positive
  expect_error(ebirdst_habitat(path, ext = e, n_predictors = 0))
  expect_error(ebirdst_habitat(path, ext = e, n_predictors = -1))
})

test_that("ebirdst_habitat pland_only", {
  predictors <- unique(habitat$predictor)
  expect_true(all(!stringr::str_detect(predictors, "(ed|sd)$")))

  # include variance metrics
  h_ed <- ebirdst_habitat(path, ext = e, pland_only = FALSE, n_predictors = 20)
  expect_length(unique(h_ed$predictor), 20)
  predictors <- unique(h_ed$predictor)
  expect_true(any(stringr::str_detect(predictors, "(ed|sd)$")))
})

test_that("plot ebirdst_habitat", {
  expect_silent({g <- plot(habitat)})
  expect_is(g, "gg")
  expect_silent({g <- plot(habitat, by_cover_class = FALSE)})
})
