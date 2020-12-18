# clean up
unlink(list.files("man", full.names = TRUE))
devtools::clean_vignettes()
pkgdown::clean_site()

# rebuild docs and install
devtools::document()
devtools::build()

# local tests and checks
devtools::test()
devtools::check(run_dont_test = TRUE)

# vignettes, readme, site
Sys.setenv(BUILD_VIGNETTES = TRUE)
devtools::build_vignettes()
rmarkdown::render("README.Rmd")
pkgdown::build_site()
Sys.unsetenv("BUILD_VIGNETTES")

# checks
devtools::check_win_devel()
devtools::check_win_release()
rhub::check_for_cran()
