library(tidyverse)
library(jsonlite)

p <- file.path("data-raw", "config.json") %>%
  read_json(simplifyVector = TRUE)

ebirdst_weeks <- tibble(week_number = seq_along(p$SRD_DATE_VEC),
                        date = as.Date(paste(p$SRD_PRED_YEAR,
                                             p$DATE_NAMES, sep = "-")),
                        week_midpoint = p$SRD_DATE_VEC)
start_end <- map(1:52, ~ c((. - 1) / 52, . / 52)) %>%
  Reduce(rbind, .) %>%
  as_tibble() %>%
  set_names(c("week_start", "week_end"))
ebirdst_weeks <- bind_cols(ebirdst_weeks, start_end)

usethis::use_data(ebirdst_weeks, overwrite = TRUE)
