library(testthat)
library(ebirdst)

sp_path <- download_data(example_data)

test_check("ebirdst")
