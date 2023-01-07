context("Data download")

skip_on_cran()
skip_if_offline()

# only run this test if an ebirdst access key is present
key <- Sys.getenv("EBIRDST_KEY")
skip_if_not(!is.na(key) && key != "" && nchar(key) > 0,
            message = "Missing ebirdst access key")

test_that("ebirdst_download() works", {
  suppressMessages({
    files <- ebirdst_download("leafly", dry_run = TRUE)
  })

  expect_true(any(grepl("config.json$", files)))
  expect_false(any(grepl("\\.db$", files)))
  expect_false(any(grepl("web_download", files)))
  expect_false(any(grepl("leafly2", files)))

  suppressMessages({
    files_all <- ebirdst_download("leafly",
                                  tifs_only = FALSE,
                                  dry_run = TRUE)
  })
  expect_true(any(grepl("config.json$", files_all)))
  expect_true(any(grepl("\\.db$", files_all)))
  expect_false(any(grepl("web_download", files_all)))
  expect_equal(length(files_all), length(files) + 2)
  expect_false(any(grepl("leafly2", files)))
})
