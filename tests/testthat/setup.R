# the contents of this file are run once for all tests

# set download directory to a temporary directory
old_dir <- Sys.getenv("EBIRDST_DATA_DIR")
temp_dir <- file.path(tempdir(), "ebirdst_temp_dir")
temp_dir <- file.path("~/Desktop/", "ebirdst_temp_dir")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(EBIRDST_DATA_DIR = temp_dir)

# download example data
path <- ebirdst_download("example_data",
                         tifs_only = FALSE,
                         show_progress = FALSE)

# cleanup the mess we made above
cleanup <- function() {
  Sys.setenv(EBIRDST_DATA_DIR = old_dir)
  unlink(temp_dir)
}
# run cleanup after tests are complete
withr::defer(cleanup(), teardown_env())

