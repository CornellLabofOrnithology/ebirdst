context("Mapping functions")
context("calc_full_extent")

test_that("ebirdst calc_full_extent", {
  library(sp)

  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82,
                    t.min = 0.425,
                    t.max = 0.475)

  abund <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                "_hr_2016_abundance_umean.tif"))
  abund <- raster_st_subset(abund, ne_extent)

  # expected RasterStack
  expect_is(calc_full_extent(abund, sp_path), "Extent")

  # RasterLayer
  expect_is(calc_full_extent(abund[[2]], sp_path), "Extent")

  # projected
  mollweide <- CRS("+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84")
  abund_prj <- suppressWarnings(raster::projectRaster(abund[[2]],
                                                      crs = mollweide))
  expect_is(calc_full_extent(abund_prj, sp_path), "Extent")

  # not a raster object
  ext_poly <- methods::as(raster::extent(abund), "SpatialPolygons")
  raster::crs(ext_poly) <- raster::crs(abund)
  expect_error(calc_full_extent(ext_poly, sp_path), "Input must be a Raster")
})

context("calc_bins")

test_that("ebirdst calc_bins", {
  library(sp)

  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82,
                    t.min = 0.425,
                    t.max = 0.475)

  abund <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                "_hr_2016_abundance_umean.tif"))
  abund <- raster_st_subset(abund, ne_extent)

  # expect a list greater than 1 for RasterStack
  expect_gt(length(calc_bins(abund)), 1)

  # expect a list greater than 1 for RasterLayer
  expect_gt(length(calc_bins(abund[[2]])), 1)

  # projected
  abund_prj <- suppressWarnings(raster::projectRaster(abund[[2]],
                                                      crs = "+init=epsg:4326"))
  expect_gt(length(calc_bins(abund_prj)), 1)

  # not a raster object
  ext_poly <- methods::as(raster::extent(abund), "SpatialPolygons")
  raster::crs(ext_poly) <- raster::crs(abund)
  expect_error(calc_bins(ext_poly), "Input must be a Raster")

  # if all is NA
  abund_test <- abund[[2]]
  abund_test[!is.na(abund_test)] <- NA
  expect_error(calc_bins(abund_test), "must have non-NA values")
})

context("map_centroids")

test_that("ebirdst map_centroids", {
  pis <- load_pis(sp_path)
  pds <- load_pds(sp_path)

  # expected without st_extent
  expect_error(map_centroids(pis = pis, pds = pds), NA)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(map_centroids(pis = pis, pds = pds, st_extent = ne_extent), NA)

  # param variations
  expect_error(map_centroids(pis = pis,
                             pds = pds,
                             st_extent = ne_extent,
                             plot_pis = FALSE), NA)
  expect_error(map_centroids(pis = pis,
                             pds = pds,
                             st_extent = ne_extent,
                             plot_pds = FALSE), NA)
  expect_error(map_centroids(pis = pis,
                             st_extent = ne_extent,
                             plot_pds = FALSE), NA)
  expect_error(map_centroids(pds = pds,
                             st_extent = ne_extent,
                             plot_pis = FALSE), NA)

  # param mis-setting
  expect_error(map_centroids(pis = pis,
                             st_extent = ne_extent,
                             plot_pis = FALSE), 'argument "pds" is missing')

  expect_error(map_centroids(pds = pds,
                             st_extent = ne_extent,
                             plot_pds = FALSE), 'argument "pis" is missing')
  expect_error(map_centroids(pis = pis,
                             pds = pds,
                             st_extent = ne_extent,
                             plot_pis = FALSE,
                             plot_pds = FALSE), "Plotting of both PIs and PDs")
})

context("calc_effective_extent")

test_that("ebirdst calc_effective_extent", {
  pis <- load_pis(sp_path)
  pds <- load_pds(sp_path)

  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 46,
                    lon.min = -86,
                    lon.max = -82,
                    t.min = 0.425,
                    t.max = 0.475)

  # expected
  expect_is(calc_effective_extent(st_extent = ne_extent, pis = pis,
                                  path = sp_path),
            "RasterLayer")
  expect_error(calc_effective_extent(st_extent = ne_extent, pis = pis,
                                     path = sp_path),
               NA)

  # known input error checking
  expect_error(calc_effective_extent(st_extent = ne_extent, pis = pis,
                                     pds = pds, path = sp_path),
               "Unable to calculate for both PIs and PDs")
  expect_error(calc_effective_extent(st_extent = ne_extent, path = sp_path),
               "Both PIs and PDs are NA")

  # no temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -86,
                    lon.max = -83)
  expect_warning(calc_effective_extent(st_extent = ne_extent, pis = pis,
                                       path = sp_path),
                 "Without temporal limits")
})
