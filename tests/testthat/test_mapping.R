context("Mapping functions")
skip_on_cran()
sp_path <- download_data("example_data", tifs_only = FALSE)
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.5, 0.6))

context("calc_full_extent")

test_that("ebirdst calc_full_extent", {
  abund <- load_raster("abundance_umean", sp_path)
  abund <- ebirdst_subset(abund, lp_extent)

  # expected RasterStack
  expect_is(calc_full_extent(abund), "Extent")

  # RasterLayer
  expect_is(calc_full_extent(abund[[2]]), "Extent")

  # projected
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
  abund_prj <- suppressWarnings(raster::projectRaster(abund[[2]],
                                                      crs = mollweide))
  expect_is(calc_full_extent(abund_prj), "Extent")

  # not a raster object
  bb <- sf::st_bbox(abund)
  expect_error(calc_full_extent(bb))
})

context("calc_bins")

test_that("ebirdst calc_bins", {
  abund <- load_raster("abundance_umean", sp_path)
  abund <- ebirdst_subset(abund, lp_extent)

  # expect a list greater than 1 for RasterStack
  expect_length(calc_bins(abund), 2)
  expect_named(calc_bins(abund), c("bins", "power"))

  # expect a list greater than 1 for RasterLayer
  expect_length(calc_bins(abund[[2]]), 2)
  expect_named(calc_bins(abund[[2]]), c("bins", "power"))

  # projected
  mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
  abund_prj <- suppressWarnings(raster::projectRaster(abund[[2]],
                                                      crs = mollweide,
                                                      method = "ngb"))
  expect_length(calc_bins(abund_prj), 2)
  expect_named(calc_bins(abund_prj), c("bins", "power"))

  # not a raster object
  bb <- sf::st_bbox(abund)
  expect_error(calc_bins(bb))

  # if all is NA
  abund_test <- abund[[2]]
  suppressWarnings(abund_test[] <- NA)
  expect_error(calc_bins(abund_test))
})

context("map_centroids")

test_that("ebirdst map_centroids", {
  # expected without st_extent
  expect_error(map_centroids(sp_path))

  # expected with ebirdst_extent
  expect_error(map_centroids(sp_path, lp_extent, plot_pis = FALSE,
                             plot_pds = FALSE))

  # returns nothing
  expect_null(map_centroids(sp_path, lp_extent))
})

context("calc_effective_extent")

test_that("ebirdst calc_effective_extent", {
  # expected
  expect_is(calc_effective_extent(sp_path, lp_extent, plot = FALSE),
            "RasterLayer")

  # known input error checking
  expect_error(calc_effective_extent(sp_path, lp_extent, "both"))
  expect_error(calc_effective_extent("/bad/path/", lp_extent, "both"))

  # no temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_warning(calc_effective_extent(sp_path, nt))
})
