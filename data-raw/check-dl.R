library(ebirdst)
library(tidyverse)

dl_and_check <- function(s, d = "~/Desktop/ebirdst/") {
  if (any(str_detect(list.files(d, "rds$"), s))) {}
  dl_path <- ebirdst_download(s, path = d)
  f_cfg <- list.files(dl_path, "rds$", recursive = TRUE, full.names = TRUE)
  l <- readRDS(f_cfg)
  l$tif_list <- list.files(dl_path, "seasonal.*tif$", recursive = TRUE) %>%
    basename()
  file.copy(f_cfg, file.path(d, basename(f_cfg)), overwrite = TRUE)

  # checks
  checks <- tibble(species_code = s, seasonal_check = TRUE, fac_check = TRUE)
  if (length(l$tif_list) > 5) {
    checks$seasonal_check <- FALSE
    message(str_glue("{s}: too many seasonal tiffs"))
  } else if (any(str_detect(l$tif_list, "breeding")) &&
      any(str_detect(l$tif_list, "year"))) {
    checks$seasonal_check <- FALSE
    message(str_glue("{s}: seasonal and year round tiffs present"))
  }
  needs <- c("CUS_PRJ", "FA_EXTENT", "SINU_EXTENT", "ABUND_BINS")
  if (!all(needs %in% names(l))) {
    checks$fac_check <- FALSE
    message(str_glue("{s}: missing fac map parameters"))
  }
  unlink(dl_path, recursive = TRUE)

  # good!
  if (checks$fac_check && checks$seasonal_check) {
    message(message(str_glue(">>>> {s}: good")))
  }
  return(checks)
}
checks <- map_dfr(ebirdst_runs$species_code, dl_and_check)
