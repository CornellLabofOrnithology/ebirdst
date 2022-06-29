context("Abundance palette")

test_that("abundance_palette", {
  expect_is(abundance_palette(n = 10), "character")
  expect_length(abundance_palette(n = 10), 10)
  expect_length(abundance_palette(n = 10, season = "breeding"), 10)
  expect_match(abundance_palette(n = 10, season = "nonbreeding"),
               "#[0-9A-F]{6}")
})
