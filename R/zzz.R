.onAttach <- function(libname, pkgname) {
  m <- paste("Please cite the eBird Status & Trends data using:",
             "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S. Ligocki, B. Petersen, ",
             "C. Wood, I. Davies, B. Sullivan, M. Iliff, S. Kelling. 2019. eBird Status and Trends, ",
             "Version: November 2019. Cornell Lab of Ornithology, Ithaca, New York.",
             "https://doi.org/10.2173/ebirdst.2019",
             sep = " \n  ")
  packageStartupMessage(m)
}
