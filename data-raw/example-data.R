library(readr)
repo <- "https://raw.githubusercontent.com/ebird/ebirdst_example-data/main/"
read_lines(file.path(repo, "example-data/2021/file-list.txt")) %>%
  write_lines("inst/extdata/example-data_file-list.txt")
