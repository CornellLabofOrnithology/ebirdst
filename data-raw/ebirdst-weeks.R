library(tidyverse)

l <- file.path("data-raw", "config.rds") %>%
  readRDS()
ebirdst_weeks <- tibble(week_number = seq_along(l$SRD_DATE_VEC),
                        date = as.Date(paste(l$SRD_PRED_YEAR,
                                             l$DATE_NAMES, sep = "-")),
                        week_midpoint = l$SRD_DATE_VEC)
start_end <- map(1:52, ~ c((. - 1) / 52, . / 52)) %>%
  Reduce(rbind, .) %>%
  as_tibble() %>%
  set_names(c("week_start", "week_end"))
ebirdst_weeks <- bind_cols(ebirdst_weeks, start_end)

usethis::use_data(ebirdst_weeks, overwrite = TRUE)
