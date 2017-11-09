context("Loading functions")

test_that("stemhelper stack_stem default", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  abund <- stack_stem(sp_path, variable = "abundance_umean")
  expect_is(abund, "RasterStack")

  expect_error(stack_stem(sp_path, variable = "misspelled"),
               paste("Selected variable is not one of the following: ",
                     "abundance_ensemble_support, abundance_lower, ",
                     "abundance_upper, abundance_umean, occurrence_umean.",
                     sep = ""))

  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(stack_stem(sp_path, variable = "abundance_umean"))
})

test_that("stemhelper stack_stem st_extent", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  abund <- stack_stem(sp_path,
                      variable = "abundance_umean",
                      st_extent = ne_extent)
  expect_is(abund, "RasterStack")
  expect_gt(raster::ncell(abund), 1)

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_is(stack_stem(sp_path,
                       variable = "abundance_umean",
                       st_extent = ne_extent),
            "RasterStack")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(stack_stem(sp_path,
                       variable = "abundance_umean",
                       st_extent = ne_extent))

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80)
  expect_error(stack_stem(sp_path,
                          variable = "abundance_umean",
                          st_extent = ne_extent))

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(stack_stem(sp_path,
                          variable = "abundance_umean",
                          st_extent = ne_extent))

  # extent with shapefile
  library(sp)
  e_extent <- raster::extent(-80, -70, 40, 47)
  e_polys <- methods::as(e_extent, "SpatialPolygons")
  raster::crs(e_polys) <- sp::CRS("+init=epsg:4326")

  ne_extent <- list(type = "polygon",
                    polygon = e_polys,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_is(stack_stem(sp_path,
                       variable = "abundance_umean",
                       st_extent = ne_extent), "RasterStack")
  expect_gt(raster::ncell(stack_stem(sp_path,
                                     variable = "abundance_umean",
                                     st_extent = ne_extent)), 1)

})

test_that("stemhelper stack_stem w/ use_analysis_extent", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  abund <- stack_stem(sp_path,
                      variable = "abundance_umean",
                      use_analysis_extent = FALSE)
  expect_is(abund, "RasterStack")
  expect_gt(raster::ncell(abund), 1)
})
