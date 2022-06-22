context("PPM functions")

skip_on_cran()

e <- ebirdst_extent(c(xmin = -90, xmax = -82, ymin = 41, ymax = 48),
                    t = c(0.25, 0.75))
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
