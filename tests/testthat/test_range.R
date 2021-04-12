test_that("make_range_polygon()", {
  r <- raster(r)

  # empty raster
  expect_null(make_range_polygon(r, layer = 1))

  # single layer, single block
  r[1:4] <- 145
  p <- make_range_polygon(r, layer = 1)
  expect_s3_class(p, "sf")
  expect_equal(nrow(p), 1)

  # single layer, two blocks of cells
  r[1:4] <- 145
  r[95:100] <- 23
  p <- make_range_polygon(r, layer = 1)
  expect_s3_class(p, "sf")
  expect_equal(nrow(p), 1)

  # all zeros
  r <- raster(r)
  values(r) <- 0
  p <- suppressWarnings(make_range_polygon(r, layer = 1, type = "range"))
  expect_null(p)
  p <- make_range_polygon(r, layer = 1, type = "prediction_area")
  expect_s3_class(p, "sf")
  expect_equal(nrow(p), 1)
  e <- extent(r) %>%
    extent_to_sf(crs = projection(r)) %>%
    st_bbox()
  expect_equal(st_bbox(p), e)

  # multiple layers, reference by names
  r1 <- raster(r)
  values(r1) <- 0
  r2 <- raster(r)
  values(r2) <- 56
  s <- stack(r1, r2) %>%
    setNames(c("breeding", "nonbreeding"))
  p <- suppressWarnings(make_range_polygon(s, layer = 1, type = "range"))
  expect_null(p)

  p <- make_range_polygon(s, layer = "nonbreeding", type = "range")
  expect_s3_class(p, "sf")
  expect_equal(nrow(p), 1)
  expect_equal(st_bbox(p), e)

  p <- make_range_polygon(s, layer = "breeding", type = "prediction_area")
  expect_s3_class(p, "sf")
  expect_equal(nrow(p), 1)
  expect_equal(st_bbox(p), e)

  # error when an invalid layer is supplied
  expect_error(make_range_polygon(s, layer = 0))
  expect_error(make_range_polygon(s, layer = 3))
  expect_error(make_range_polygon(s, layer = "migration"))
})
