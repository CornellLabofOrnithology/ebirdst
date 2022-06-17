.onAttach <- function(libname, pkgname) {
  m <- paste("Please cite the eBird Status & Trends data using:",
             "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S. Ligocki, W. Hochachka,",
             "C. Wood, I. Davies, M. Iliff, L. Seitz. %s. eBird Status and Trends,",
             "Data Version: %s; Released: %s Cornell Lab of Ornithology, Ithaca, New York.",
             "https://doi.org/10.2173/ebirdst.%s",
             "\nNOTE: after November 2022, the %s version of the data will no longer be available.",
             "Those wishing to continue using Status & Trends data will be required to transition",
             "to using a newer version.",
             sep = " \n  ")
  v <- ebirdst_version()
  dv <- v["data_version"]
  rv <- v["release_year"]
  packageStartupMessage(sprintf(m, rv, dv, rv, dv, dv))
}
