library(testthat)
library(ebirdst)

species <- "yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83"
sp_path <- download_data("~/tmp/", species = species, example_data = TRUE)

test_check("ebirdst")
