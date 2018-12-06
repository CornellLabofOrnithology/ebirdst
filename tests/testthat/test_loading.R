context("Loading functions")
context("raster_st_subset")

test_that("ebirdst raster_st_subset", {
  abunds <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                 "_hr_2016_abundance_umean.tif"))
  expect_is(abunds, "RasterStack")

  ne_extent <- list(type = "rectangle",
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82,
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
                    lat.min = 42,
                    lat.max = 47,
                    lon.min = -88,
                    lon.max = -82)
  expect_warning(raster_st_subset(abunds, ne_extent))

  abund_sub <- raster_st_subset(abunds, ne_extent)

  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 52)
  expect_equal(class(abund_sub)[1], "RasterBrick")

  # extent with shapefile
  library(sp)
  e_extent <- raster::extent(-88, -82, 40, 47)
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

context("label_raster_stack and parse_raster_dates")

test_that("ebirdst label_raster_stack", {
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

  # expected
  date_vector <- parse_raster_dates(abunds_labeled)
  expect_equal(class(date_vector), "Date")

  # error
  expect_error(parse_raster_dates(abunds),
               "raster names not in correct format. Please call label_raster_stack first.")
})

context("load_pis")

test_that("ebirdst load_pis", {
  expect_is(load_pis(sp_path), "data.frame")
  expect_gt(nrow(load_pis(sp_path)), 0)

  # broken path
  expect_error(load_pis('~/some/messed/up/path/that/does/not/exist'),
               "RData file does not exist")
})

context("load_pds")

test_that("stemhelpe load_pds", {
  expect_is(load_pds(sp_path), "data.frame")
  expect_gt(nrow(load_pds(sp_path)), 0)

  # broken path
  expect_error(load_pds('~/some/messed/up/path/that/does/not/exist'),
               "RData file does not exist")
})
