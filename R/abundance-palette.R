#' eBird Status and Trends color palettes for mapping
#'
#' Generate the color palettes used for the eBird Status and Trends relative
#' abundance maps.
#'
#' @param n integer; the number of colors to be in the palette.
#' @param season character; the season to generate colors for or "weekly" to
#'   get the color palette used in the weekly abundance animations.
#'
#' @return A character vector of hex color codes.
#' @export
#'
#' @examples
#' # breeding season color palette
#' abundance_palette(10, season = "breeding")
abundance_palette <- function(n, season = c("weekly",
                                            "breeding", "nonbreeding",
                                            "migration",
                                            "prebreeding_migration",
                                            "postbreeding_migration",
                                            "year_round")) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  season <- match.arg(season)

  # set base color by season
  col_zero <- "#e6e6e6"
  if (season == "weekly") {
    plsm <- rev(viridisLite::plasma(n - 1, end = 0.9))
    plsm <- stringr::str_remove(plsm, "FF$")
    gry <- grDevices::colorRampPalette(c(col_zero, plsm[1]))
    return(c(gry(4)[2], plsm))
  } else if (season == "breeding") {
    base_col <- "#cc503e"
  } else if (season == "nonbreeding") {
    base_col <- "#1d6996"
  } else if (season %in% c("migration", "postbreeding_migration")) {
    base_col <- "#edad08"
  } else if (season == "prebreeding_migration") {
    base_col <- "#73af48"
  } else if (season == "year_round") {
    base_col <- "#6f4070"
  } else {
    stop("Invalid season.")
  }

  # seasonal palettes
  gry <- grDevices::colorRampPalette(c(col_zero, base_col))
  mid <- grDevices::colorRampPalette(c(gry(5)[2], base_col))
  black <- grDevices::colorRampPalette(c(base_col, "#000000"))
  pal <- grDevices::colorRampPalette(c(gry(5)[2], mid(9)[5], base_col,
                                       black(5)[2]))
  return(pal(n))
}
