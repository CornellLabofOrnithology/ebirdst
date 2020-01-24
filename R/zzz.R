.onAttach <- function(libname, pkgname) {
  m <- paste("Please cite the eBird Status & Trends data using:",
             "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S. Ligocki, B. Petersen, ",
             "C. Wood, I. Davies, B. Sullivan, M. Iliff, S. Kelling. 2019. eBird Status and Trends, ",
             "Version: November 2019. https://ebird.org/science/status-and-trends. Cornell Lab of ",
             "Ornithology, Ithaca, New York. https://doi.org/...",
             sep = " \n  ")
  packageStartupMessage(m)
}
