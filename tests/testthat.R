library(testthat)
library(ebirdst)

dl_dir <- tempdir()
sp_path <- download_data(example_data, path = dl_dir)

test_check("ebirdst")
