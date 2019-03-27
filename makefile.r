unlink(file.path(rappdirs::user_data_dir("ebirdst"),
          "yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83"),
       recursive = TRUE)
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
file.copy(list.files(".", "README.*png$", full.names = TRUE), "docs/")
Sys.unsetenv("BUILD_VIGNETTES")

chk <- rhub::check_for_cran(platforms = c("windows-x86_64-devel",
                                          "debian-gcc-devel",
                                          "fedora-clang-devel"))
chk
