.onAttach <- function(libname, pkgname) {
  v <- ebirdst_version()
  vy <- v[["version_year"]]
  ry <- v[["release_year"]]
  m <- stringr::str_glue(
    "Please cite the eBird Status & Trends data using: ",
    "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S. Ligocki,",
    "W. Hochachka, L. Jaromczyk, C. Wood, I. Davies, M. Iliff, L. Seitz. {ry}.",
    "eBird Status and Trends, Data Version: {vy}; Released: {ry}. Cornell Lab of",
    "Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.{vy}",
    .sep = "\n  ", .trim = FALSE
  )
  packageStartupMessage(as.character(m))
}
