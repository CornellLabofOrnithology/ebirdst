context("PPM functions")

# compute ppms
# plot binary by time
# plot all ppms

test_that("stemhelper compute_ppms", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  # expected
  expect_error(compute_ppms(sp_path), NA)

  cppms <- compute_ppms(sp_path)

  expect_equal(length(cppms), 3)
  expect_is(cppms, "list")
  expect_is(cppms$binary_stats, "data.frame")
  expect_is(cppms$occ_stats, "data.frame")
  expect_is(cppms$count_stats, "data.frame")
  expect_is(cppms$binary_stats$mc, "integer")
  expect_is(cppms$occ_stats$mc, "integer")
  expect_is(cppms$count_stats$mc, "integer")
  expect_equal(length(cppms$binary_stats$mc), 25)
  expect_equal(length(cppms$occ_stats$mc), 25)
  expect_equal(length(cppms$count_stats$mc), 25)

  # with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent), NA)

  cppms <- compute_ppms(sp_path, st_extent = ne_extent)

  expect_equal(length(cppms), 3)
  expect_is(cppms, "list")
  expect_is(cppms$binary_stats, "data.frame")
  expect_is(cppms$occ_stats, "data.frame")
  expect_is(cppms$count_stats, "data.frame")
  expect_is(cppms$binary_stats$mc, "integer")
  expect_is(cppms$occ_stats$mc, "integer")
  expect_is(cppms$count_stats$mc, "integer")
  expect_equal(length(cppms$binary_stats$mc), 25)
  expect_equal(length(cppms$occ_stats$mc), 25)
  expect_equal(length(cppms$count_stats$mc), 25)

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent), NA)

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(compute_ppms(sp_path, st_extent = ne_extent),
               "st_extent argument must be a list")
})
