context("PPM functions")

skip_on_cran()
skip_on_appveyor()
skip_on_travis()

path <- ebirdst_download("example_data", tifs_only = FALSE)
e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43),
                    t = c(0.5, 0.6))
ppm <- ebirdst_ppms(path, ext = e)
wk <- ebirdst_ppms_ts(path, ext = e, summarize_by = "weeks")
mt <- ebirdst_ppms_ts(path, ext = e, summarize_by = "months")

test_that("ebirdst_ppms", {
  # expected
  expect_equal(length(ppm), 3)
  expect_is(ppm, "ebirdst_ppms")
  expect_named(ppm, c("binary_ppms", "occ_ppms", "abd_ppms"))
  expect_is(ppm$binary_ppms, "data.frame")
  expect_is(ppm$occ_ppms, "data.frame")
  expect_is(ppm$abd_ppms, "data.frame")
  expect_is(ppm$binary_ppms$mc_iteration, "integer")
  expect_is(ppm$occ_ppms$mc_iteration, "integer")
  expect_is(ppm$abd_ppms$mc_iteration, "integer")
  expect_length(ppm$binary_ppms$mc_iteration, 25)
  expect_length(ppm$occ_ppms$mc_iteration, 25)
  expect_length(ppm$abd_ppms$mc_iteration, 25)

  # invalid inputs
  expect_error(ebirdst_ppms("/invalid/path", ext = e))
})

test_that("ebirdst_ppms_ts", {
  # weeks
  expect_equal(length(wk), 3)
  expect_is(wk, "ebirdst_ppms_ts")
  expect_named(wk, c("binary_ppms", "occ_ppms", "abd_ppms"))
  for (t in c("binary_ppms", "occ_ppms", "abd_ppms")) {
    expect_is(wk[[t]], "data.frame")
    expect_is(wk[[t]][["week"]], "Date")
    expect_true(all(ebirdst::ebirdst_weeks$date %in% wk[[t]][["week"]]))
  }

  # months
  expect_equal(length(mt), 3)
  expect_is(mt, "ebirdst_ppms_ts")
  expect_named(mt, c("binary_ppms", "occ_ppms", "abd_ppms"))
  for (t in c("binary_ppms", "occ_ppms", "abd_ppms")) {
    expect_is(mt[[t]], "data.frame")
    expect_is(mt[[t]][["month"]], "factor")
    expect_true(all(month.abb %in% mt[[t]][["month"]]))
  }

  # invalid inputs
  expect_error(ebirdst_ppms_ts("/invalid/path", ext = e))
  expect_error(ebirdst_ppms_ts(path, ext = e, summarize_by = "days"))
  expect_error(ebirdst_ppms_ts(path, ext = data.frame()))
})

test_that("plot ebirdst_ppms", {
  expect_silent(plot(ppm))
})

test_that("plot ebirdst_ppms_ts", {
  # weeks
  expect_silent(plot(wk))
  # binary, kappa
  expect_silent(plot(wk, type = "binary", metric = "kappa"))
  # occurrence, sensitivity
  expect_silent(plot(wk, type = "occurrence", metric = "sensitivity"))
  # abundance, poisson deviance
  expect_silent(plot(wk, type = "abundance", metric = "poisson_dev_abd"))

  # months
  expect_silent(plot(mt))
  # binary, kappa
  expect_silent(plot(mt, type = "binary", metric = "kappa"))
  # occurrence, sensitivity
  expect_silent(plot(mt, type = "occurrence", metric = "sensitivity"))
  # abundance, poisson deviance
  expect_silent(plot(mt, type = "abundance", metric = "poisson_dev_abd"))

  # invalid inputs
  expect_error(plot(wk, type = "invalid", metric = "kappa"))
  expect_error(plot(mt, type = "invalid", metric = "kappa"))
  expect_error(plot(wk, type = "binary", metric = "poisson_dev_abd"))
  expect_error(plot(mt, type = "abundance", metric = "kappa"))
})
