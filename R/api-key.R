#' Store the eBird Status and Trends API key
#'
#' Accessing eBird Status and Trends data requires an API key, which can be
#' obtained by visiting https://ebird.org/st/request. This key must be
#' stored as the environment variable `EBIRDST_KEY` in order for
#' [ebirdst_download()] to use it. The easiest approach is to store the key in
#' your `.Renviron` file so it can always be accessed in your R sessions. Use
#' this function to set `EBIRDST_KEY` in your `.Renviron` file provided that it
#' is located in the standard location in your home directory. It is also
#' possible to manually edit the `.Renviron` file. **The API key is specific to
#' you and should never be shared or made publicly accessible.**
#'
#' @param key character; API key obtained by filling out the form at
#'   https://ebird.org/st/request.
#' @param overwrite logical; should the existing `EBIRDST_KEY` be overwritten if
#'   it has already been set in .Renviron.
#'
#' @return Edits .Renviron, then returns the path to this file invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # save the api key, replace XXXXXX with your actual key
#' set_ebirdst_api_key("XXXXXX")
#' }
set_ebirdst_api_key <- function(key, overwrite = FALSE) {
  stopifnot(is.character(key), length(key) == 1, nchar(key) > 0)
  stopifnot(is.logical(overwrite), length(overwrite) == 1)

  # find .Renviron
  renv_path <- path.expand(file.path("~", ".Renviron"))
  if (!file.exists(renv_path)) {
    file.create(renv_path)
  }
  renv_lines <- readLines(renv_path)

  key_line <- paste0("EBIRDST_KEY='", key, "'")

  # look for existing entry, remove if overwrite = TRUE
  renv_exists <- grepl("^EBIRDST_KEY[[:space:]]*=.*", renv_lines)
  if (sum(renv_exists) > 1) {
    stop("Multiple entries for EBIRDST_KEY appear in .Renviron")
  } else if (any(renv_exists)) {
    if (isTRUE(overwrite)) {
      # replace existing
      renv_lines[which(renv_exists)] <- key_line
    } else {
      stop("EBIRDST_KEY already set, use overwrite = TRUE to overwite.")
    }
  } else {
    renv_lines <- c(renv_lines, key_line)
  }
  # set key in .Renviron
  writeLines(renv_lines, renv_path)
  invisible(renv_path)
}


get_ebirdst_api_key <- function() {
  key <- Sys.getenv("EBIRDST_KEY")
  if (is.na(key) || key == "" || nchar(key) == 0) {
    stop("eBird Status and Trends API key not found in .Renviron. ",
         "Considering using set_ebirdst_api_key() to store the key.")
  }
  invisible(key)
}
