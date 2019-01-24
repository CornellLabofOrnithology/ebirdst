unlink(list.files("man", full.names = TRUE))
devtools::clean_vignettes()
pkgdown::clean_site()

devtools::document()

devtools::install_local(force = TRUE)

devtools::test()
devtools::check(run_dont_test = TRUE)

Sys.setenv(BUILD_VIGNETTES = TRUE)
devtools::build_vignettes()
rmarkdown::render("README.Rmd")
pkgdown::build_site()
file.copy(list.files(".", "README.*png", full.names = TRUE), "docs/")
Sys.unsetenv("BUILD_VIGNETTES")
