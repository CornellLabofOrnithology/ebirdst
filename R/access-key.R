#' Store the eBird Status and Trends access key
#'
#' Accessing eBird Status and Trends data requires an access key, which can be
#' obtained by visiting https://ebird.org/st/request. This key must be
#' stored as the environment variable `EBIRDST_KEY` in order for
#' [ebirdst_download()] to use it. The easiest approach is to store the key in
#' your `.Renviron` file so it can always be accessed in your R sessions. Use
#' this function to set `EBIRDST_KEY` in your `.Renviron` file provided that it
#' is located in the standard location in your home directory. It is also
#' possible to manually edit the `.Renviron` file. **The access key is specific
#' to you and should never be shared or made publicly accessible.**
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
#' set_ebirdst_access_key("XXXXXX")
#' }
set_ebirdst_access_key <- function(key, overwrite = FALSE) {
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
    # set key in .Renviron
    writeLines(renv_lines, renv_path)
  } else {
    write(paste0("\n", key_line, "\n"), renv_path, append = TRUE)
  }
  message("eBird Status and Trends access key stored in: ", renv_path,
          "\nYou must RESTART R to load the saved access key.")
  invisible(renv_path)
}


get_ebirdst_access_key <- function() {
  key <- Sys.getenv("EBIRDST_KEY")
  if (is.na(key) || key == "" || nchar(key) == 0) {
    message("An access key is required to download eBird Status and Trends ",
            "data\n1. Get a key by filling out the request form at ",
            "https://ebird.org/st/request\n",
            "2. Save the key using set_ebirdst_access_key()\n",
            "3. Restart R to load the key")
    stop("Valid eBird Status and Trends access key not found. ",
         "Note that keys expire after 1 month, you may need a new key.")
  }
  invisible(key)
}
