context("Loading functions")
context("stack_stem")

test_that("stemhelper stack_stem default", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)

  abund <- stack_stem(sp_path, variable = "abundance_umean", st_extent = ne_extent)
  expect_is(abund, "RasterStack")

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.435)

  abund <- stack_stem(sp_path, variable = "abundance_umean",
                      st_extent = ne_extent)
  expect_is(abund, "RasterLayer")

  expect_error(stack_stem(sp_path, variable = "misspelled"),
               paste("Selected variable is not one of the following: ",
                     "abundance_lower, ",
                     "abundance_upper, abundance_umean, occurrence_umean.",
                     sep = ""))

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(stack_stem(sp_path, variable = "abundance_umean"),
               "RData file does not exist")
})

test_that("stemhelper stack_stem st_extent", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.435)
  abund <- stack_stem(sp_path, variable = "abundance_umean",
                      st_extent = ne_extent)
  expect_is(abund, "RasterLayer")
  expect_gt(raster::ncell(abund), 1)

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_is(stack_stem(sp_path, variable = "abundance_umean",
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
  expect_error(stack_stem(sp_path, variable = "abundance_umean",
                          st_extent = ne_extent),
               "Minimum latitude is greater than maximum latitude")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80)
  expect_error(stack_stem(sp_path, variable = "abundance_umean",
                          st_extent = ne_extent),
               "Either lon.min or lon.max missing")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(stack_stem(sp_path, variable = "abundance_umean",
                          st_extent = ne_extent),
               "st_extent argument must be a list")

  # extent with shapefile
  library(sp)
  e_extent <- raster::extent(-80, -70, 40, 47)
  e_polys <- methods::as(e_extent, "SpatialPolygons")
  raster::crs(e_polys) <- sp::CRS("+init=epsg:4326")

  ne_extent <- list(type = "polygon",
                    polygon = e_polys,
                    t.min = 0.425,
                    t.max = 0.435)
  abund <- stack_stem(sp_path, variable = "abundance_umean",
                      st_extent = ne_extent)
  expect_is(abund, "RasterLayer")
  expect_gt(raster::ncell(abund), 1)

})

test_that("stemhelper stack_stem w/ use_analysis_extent", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.435)
  abund <- stack_stem(sp_path,
                      variable = "abundance_umean",
                      st_extent = ne_extent,
                      use_analysis_extent = FALSE)
  expect_is(abund, "RasterLayer")
  expect_gt(raster::ncell(abund), 1)
})

context("load_pis")

test_that("stemhelper load_pis", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  expect_is(load_pis(sp_path), "data.frame")
  expect_gt(nrow(load_pis(sp_path)), 0)

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(load_pis(sp_path), "RData file does not exist")
})

context("load_pds")

test_that("stemhelpe load_pds", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  expect_is(load_pds(sp_path), "data.frame")
  expect_gt(nrow(load_pds(sp_path)), 0)

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(load_pds(sp_path), "RData file does not exist")
})
