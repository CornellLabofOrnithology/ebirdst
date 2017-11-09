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
