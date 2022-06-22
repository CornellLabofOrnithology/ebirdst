context("Loading ranges")

skip_on_cran()

test_that("load_ranges", {
  # raw
  range_raw <- load_ranges(path, resolution = "lr", smoothed = FALSE)
  expect_s3_class(range_raw, "sf")
  expect_equal(nrow(range_raw), 4)
  expect_named(range_raw,
               c("species_code", "scientific_name", "common_name",
                 "prediction_year", "type", "season", "start_date", "end_date",
                 "geom"))

  # smoothed
  range_smooth <- load_ranges(path, resolution = "lr")
  expect_s3_class(range_smooth, "sf")
  expect_equal(nrow(range_smooth), 4)
  expect_named(range_smooth,
               c("species_code", "scientific_name", "common_name",
                 "prediction_year", "type", "season", "start_date", "end_date",
                 "geom"))
})

test_that("load_ranges error", {
  expect_error(load_ranges(path, "error"))
  expect_error(load_ranges("/invalid/path/"))
})
