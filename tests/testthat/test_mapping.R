context("Mapping functions")

test_that("stemhelper calc_full_extent", {
  library(sp)

  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  abund <- stack_stem(sp_path, variable = "abundance_umean")

  # expected RasterStack
  expect_is(calc_full_extent(abund), "Extent")

  # RasterLayer
  expect_is(calc_full_extent(abund[[26]]), "Extent")

  # projected
  abund_prj <- raster::projectRaster(abund[[26]], crs = "+init=epsg:4326")
  expect_is(abund_prj, "Extent")

  # not a raster object
  ext_poly <- methods::as(raster::extent(abund), "SpatialPolygons")
  raster::crs(ext_poly) <- raster::crs(abund)
  expect_error(calc_full_extent(ext_poly), "Input must be a Raster")
})

test_that("stemhelper calc_bins", {
  library(sp)

  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  abund <- stack_stem(sp_path, variable = "abundance_umean")

  # expected RasterStack
  expect_gt(length(calc_bins(abund)), 1)

  # RasterLayer
  expect_gt(length(calc_bins(abund[[26]])), 1)

  # projected
  abund_prj <- raster::projectRaster(abund[[26]], crs = "+init=epsg:4326")
  expect_gt(length(calc_bins(abund_prj)), 1)

  # not a raster object
  ext_poly <- methods::as(raster::extent(abund), "SpatialPolygons")
  raster::crs(ext_poly) <- raster::crs(abund)
  expect_error(calc_bins(ext_poly), "Input must be a Raster")

  # if all is NA
  abund_test <- abund[[26]]
  abund_test[!is.na(abund_test)] <- NA
  expect_error(calc_bins(abund_test), "must have non-NA values")
})

test_that("stemhelper combine_layers", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  abund <- stack_stem(sp_path, variable = "abundance_umean")

  # expected
  expect_is(combine_layers(abund, sp_path, 26), "RasterLayer")
  expect_false(is.na(raster::maxValue(combine_layers(abund, sp_path, 26))))

  # not a RasterStack
  expect_error(combine_layers(abund[[26]], sp_path, 26),
               "stack must be of type RasterStack")

  # week is not integer
  # week outside of 1 to 52
  expect_error(combine_layers(abund, sp_path, 4.1), "week must be an integer")
  expect_error(combine_layers(abund, sp_path, 0), "week must be an integer")
  expect_error(combine_layers(abund, sp_path, 101), "week must be an integer")

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(combine_layers(abund, sp_path, 26), "RData file does not exist")
})

test_that("stemhelper map_centroids", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pis <- load_pis(sp_path)
  pds <- load_pds(sp_path)

  # expected without st_extent
  expect_error(map_centroids(pis = pis, pds = pds), NA)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
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

test_that("stemhelper calc_effective_extent", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pis <- load_pis(sp_path)
  pds <- load_pds(sp_path)

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)

  calc_effective_extent(st_extent = ne_extent, pis = pis)

  # expected
  expect_is(calc_effective_extent(st_extent = ne_extent, pis = pis),
            "RasterLayer")
  expect_is(calc_effective_extent(st_extent = ne_extent, pds = pds),
            "RasterLayer")
  expect_error(calc_effective_extent(st_extent = ne_extent, pds = pds), NA)
  expect_error(calc_effective_extent(st_extent = ne_extent, pis = pis), NA)

  # known input error checking
  expect_error(calc_effective_extent(st_extent = ne_extent,
                                     pis = pis,
                                     pds = pds),
               "Unable to calculate for both PIs and PDs")
  expect_error(calc_effective_extent(st_extent = ne_extent),
               "Both PIs and PDs are NA")

  # no temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_warning(calc_effective_extent(st_extent = ne_extent, pis = pis),
                 "Without temporal limits")
})
