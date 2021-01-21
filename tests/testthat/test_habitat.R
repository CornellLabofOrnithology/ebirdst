context("Habitat associations")

skip_on_cran()
skip_on_appveyor()
skip_on_travis()

path <- ebirdst_download("example_data", tifs_only = FALSE)
e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43))
habitat <- ebirdst_habitat(path, ext = e)

test_that("ebirdst_habitat", {
  # expected
  expect_is(habitat, "ebirdst_habitat")
  expect_is(habitat, "data.frame")
  expect_named(habitat, c("predictor", "date", "importance", "direction"))
  expect_is(habitat$predictor, "character")
  expect_is(habitat$date, "numeric")
  expect_is(habitat$importance, "numeric")
  expect_is(habitat$direction, "numeric")
  expect_length(unique(habitat$date), 52)

  # invalid inputs
  expect_error(ebirdst_habitat("/invalid/path", ext = e))
})

test_that("ebirdst_habitat extent", {
  # temporal extent ignored
  e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43),
                      t = c(0.5, 0.75))
  h2 <- ebirdst_habitat(path, ext = e)
  expect_identical(habitat, h2)

  # spatial extent required
  expect_error(ebirdst_habitat(path))
})

test_that("ebirdst_habitat pland_only", {
  predictors <- unique(habitat$predictor)
  expect_true(all(!stringr::str_detect(predictors, "(ed|sd)$")))
})

test_that("plot ebirdst_habitat", {
  expect_silent({g <- plot(habitat)})
  expect_is(g, "gg")
  expect_silent({g <- plot(habitat, group_roads = FALSE)})
})
