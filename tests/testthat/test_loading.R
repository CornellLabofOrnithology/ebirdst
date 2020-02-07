context("Loading and subsetting functions")

sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.5, 0.6))

test_that("subset raster", {
  #skip_on_cran()
  abunds <- load_raster("abundance", sp_path)
  #abunds <- raster::stack(f_dl)
  expect_is(abunds, "RasterStack")

  abund_sub <- ebirdst_subset(abunds, lp_extent)

  # expected
  # number of layers
  # check type
  # number of cells
  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 5)
  expect_is(abund_sub, "RasterBrick")

  ### variations
  # without temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_warning({abund_sub <- ebirdst_subset(abunds, nt)})

  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 52)
  expect_is(abund_sub, "RasterBrick")

  # extent with polygon
  poly_extent <- ebirdst_extent(sf::st_as_sfc(lp_extent$extent),
                                t = c(0.5, 0.6))
  abund_sub <- ebirdst_subset(abunds, poly_extent)
  expect_gt(raster::ncell(abund_sub), 1)
  expect_equal(raster::nlayers(abund_sub), 5)

  expect_is(abund_sub, "RasterBrick")
  ### error tests
  expect_error(load_raster("abundance", "/bad/path"))
  # broken path

  # reversed min max
  expect_error(ebirdst_extent(c(xmin = -83, xmax = -86, ymin = 42, ymax = 45)))
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47, ymax = 45)))

  # missing a corner
  expect_error(ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 47)))
})

context("label_raster_stack and parse_raster_dates")

test_that("ebirdst label_raster_stack", {
  skip_on_cran()
  abunds <- load_raster("abundance", sp_path)

  # expected
  expect_equal(length(grep("X2018.", names(abunds)[1])), 1)

  # error
  expect_error(label_raster_stack(abunds[[1:3]]))
  expect_error(label_raster_stack(abunds[[1]]))

  # expected
  date_vector <- parse_raster_dates(abunds)
  expect_is(date_vector, "Date")
})

context("load_pis")

test_that("ebirdst load_pis", {
  expect_is(load_pis(sp_path), "data.frame")
  expect_gt(nrow(load_pis(sp_path)), 0)

  # broken path
  expect_error(load_pis("~/some/messed/up/path/that/does/not/exist"))
})
