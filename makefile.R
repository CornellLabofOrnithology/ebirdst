# clean up
unlink(list.files("man", full.names = TRUE))

# rebuild docs and install
devtools::document()
pak::local_install()

# local tests
devtools::test()

# vignettes, readme, site
devtools::clean_vignettes()
pkgdown::clean_site()
Sys.setenv(BUILD_VIGNETTES = TRUE)
rmarkdown::render("README.Rmd")
unlink("README.html")
pkgdown::build_site()
Sys.unsetenv("BUILD_VIGNETTES")

# local checks
devtools::check()

# checks
devtools::check_win_devel()
devtools::check_win_release()
rhub::check_for_cran(platforms = "solaris-x86-patched", show_status = FALSE)
