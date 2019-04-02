# clean up
unlink(file.path(rappdirs::user_data_dir("ebirdst"),
                 "yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83"),
       recursive = TRUE)
unlink(list.files("man", full.names = TRUE))
devtools::clean_vignettes()
pkgdown::clean_site()

# rebuild docs and install
devtools::document()
devtools::install_local(force = TRUE)

# local tests and checks
devtools::test()
devtools::check(run_dont_test = TRUE)

# vignettes, readme, site
Sys.setenv(BUILD_VIGNETTES = TRUE)
devtools::build_vignettes()
rmarkdown::render("README.Rmd")
pkgdown::build_site()
file.copy(list.files(".", "README.*png$", full.names = TRUE), "docs/")
Sys.unsetenv("BUILD_VIGNETTES")

# checks
devtools::check_win_devel()
devtools::check_win_release()
rhub::check_for_cran(platforms = c("windows-x86_64-devel",
                                   "fedora-clang-devel"),
                     show_status = FALSE)
