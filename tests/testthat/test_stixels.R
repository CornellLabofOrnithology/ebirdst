context("Stixel functions")

skip_on_cran()

e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                    t = c(0.5, 0.6))
stx <- load_stixels(path)
stx_sf <- load_stixels(path, return_sf = TRUE)

test_that("load_stixels", {
  # expectations
  expect_is(stx, "data.frame")
  expect_equal(ncol(stx), 76)
  expect_gt(nrow(stx), 0)
  expect_true(all(c("stixel_id", "longitude", "latitude", "day_of_year") %in%
                    names(stx)))

  # spatial format
  expect_is(stx_sf, "sf")
  expect_is(sf::st_geometry(stx_sf), "sfc_POINT")
  expect_equal(ncol(stx_sf), 75)
  expect_gt(nrow(stx_sf), 0)
  expect_true(all(c("stixel_id", "day_of_year") %in% names(stx_sf)))

  # invalid input
  expect_error(load_stixels("/invalid/path"))
})

test_that("stixelize", {
  p <- stixelize(stx)

  # expectations
  expect_is(p, "sf")
  expect_is(sf::st_geometry(p), "sfc_POLYGON")
  expect_equal(ncol(stx) + 1, ncol(p))
  expect_equal(nrow(stx), nrow(p))
  expect_true(all(c("stixel_id", "longitude", "latitude", "day_of_year") %in%
                    names(p)))

  # invalid input
  expect_error(stixelize(dplyr::select(stx, -lat)))
})

test_that("stixelize for dateline crossing stixels", {
  stix_df <- data.frame(longitude_min = c(-170, -170),
                        longitude_max = c(170, -160),
                        latitude_min = c(0, 0),
                        latitude_max = c(10, 10))
  p <- stixelize(stix_df)

  # expectations
  expect_is(p, "sf")
  expect_equal(as.character(sf::st_geometry_type(p)),
               c("MULTIPOLYGON", "POLYGON"))
  expect_equal(ncol(stix_df) + 1, ncol(p))
  expect_equal(nrow(stix_df), nrow(p))
})
