# clean up
unlink(list.files("man", full.names = TRUE))
devtools::clean_vignettes()
pkgdown::clean_site()

# rebuild docs and install
devtools::document()
pak::pkg_install(".")

# local tests and checks
devtools::test()
devtools::check(run_dont_test = TRUE)

# vignettes, readme, site
Sys.setenv(BUILD_VIGNETTES = TRUE)
rmarkdown::render("README.Rmd")
unlink("README.html")
pkgdown::build_site()
Sys.unsetenv("BUILD_VIGNETTES")

# checks
devtools::check_win_devel()
devtools::check_win_release()
rhub::check_for_cran(platforms = "solaris-x86-patched", show_status = FALSE)
rhub::check_for_cran(platforms = "debian-gcc-release", show_status = FALSE)

