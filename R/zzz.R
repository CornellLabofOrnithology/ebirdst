.onAttach <- function(libname, pkgname) {
  v <- ebirdst_version()
  vy <- v[["version_year"]]
  ry <- v[["release_year"]]
  aed <- format(v[["access_end_date"]])
  citation <- stringr::str_glue(
    "Please cite the eBird Status & Trends data using: ",
    "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, O. Robinson, S. Ligocki,",
    "W. Hochachka, L. Jaromczyk, C. Wood, I. Davies, M. Iliff, L. Seitz. {ry}.",
    "eBird Status and Trends, Data Version: {vy}; Released: {ry}. Cornell Lab of",
    "Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.{vy}",
    .sep = "\n  ", .trim = FALSE)
  m <- stringr::str_glue(
    "{citation}",
    "",
    "This version of the package provides access to the {vy} version of the eBird",
    "Status Data Products. Access to the {vy} data will be provided until November",
    "2022, when it will be replaced by the {vy + 1} data. At that point, you will be",
    "required to update this R package and transition to using the new data.",
    .sep = "\n", .trim = FALSE
  )
  packageStartupMessage(as.character(m))
}
