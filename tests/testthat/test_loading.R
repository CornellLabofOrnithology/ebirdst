context("Loading functions")
context("raster_st_subset")

test_that("stemhelper raster_st_subset", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  abunds <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                 "_hr_2016_abundance_umean.tif"))
  expect_is(abunds, "RasterStack")

  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.5,
                    t.max = 0.6)

  abund_sub <- raster_st_subset(abunds, ne_extent)

  ### expected
  # number of layers
  # check type
  # number of cells
  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 5)
  expect_equal(class(abund_sub)[1], "RasterBrick")

  ### variations
  # without temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_warning(raster_st_subset(abunds, ne_extent))

  abund_sub <- raster_st_subset(abunds, ne_extent)

  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 52)
  expect_equal(class(abund_sub)[1], "RasterBrick")

  # extent with shapefile
  library(sp)
  e_extent <- raster::extent(-80, -70, 40, 47)
  e_polys <- methods::as(e_extent, "SpatialPolygons")
  raster::crs(e_polys) <- sp::CRS("+init=epsg:4326")

  ne_extent <- list(type = "polygon",
                    polygon = e_polys,
                    t.min = 0.5,
                    t.max = 0.6)
  abund_sub <- raster_st_subset(abunds, ne_extent)
  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 5)
  expect_equal(class(abund_sub)[1], "RasterBrick")

  ### error tests
  # broken path
  expect_error(raster::stack(paste0(sp_path, "/results/tifsWRONG/", species,
                                    "_hr_2016_abundance_umean.tif")))
  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(raster_st_subset(abunds, ne_extent),
               "Minimum latitude is greater than maximum latitude")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80)
  expect_error(raster_st_subset(abunds, ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(raster_st_subset(abunds, ne_extent),
               "st_extent argument must be a list")
})

context("label_raster_stack")

test_that("stemhelper label_raster_stack", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-SP_TEST-20180724-7ff34421"
  sp_path <- paste(root_path, species, sep = "")

  abunds <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                 "_hr_2016_abundance_umean.tif"))

  # expected
  abunds_labeled <- label_raster_stack(abunds)
  expect_equal(length(grep("X2016.", names(abunds_labeled)[1])), 1)

  # error
  expect_error(label_raster_stack(abunds[[1:3]]),
               "The raster_data object must be full stack or brick of 52")
  expect_error(label_raster_stack(abunds[[1]]),
               "The raster_data object must be either RasterStack or RasterBrick.")
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
