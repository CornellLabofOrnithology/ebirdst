library(tidyverse)
library(auk)
library(lubridate)

runs <- read_csv("data-raw/seasons.csv") %>%
  rename_all(tolower) %>%
  rename(resident_quality = full_year_quality) %>%
  filter(!is.na(run_name))

# remove species for which nothing passed review
seasons <- c("nonbreeding", "prebreeding_migration",
             "breeding", "postbreeding_migration",
             "resident")
passes <- rowSums(as.matrix(runs[, paste0(seasons, "_quality")]) > 0)
runs <- runs[!is.na(passes) & passes > 0, ]

# correctly na season dates
for (s in seasons) {
  s_fail <- runs[[paste0(s, "_quality")]] == 0 |
    is.na(runs[[paste0(s, "_quality")]])
  runs[[paste0(s, "_start")]][s_fail] <- NA_character_
  runs[[paste0(s, "_end")]][s_fail] <- NA_character_
  if (s == "resident") {
    runs[[paste0(s, "_start")]][!runs$resident] <- NA_character_
    runs[[paste0(s, "_end")]][!runs$resident] <- NA_character_
  } else {
    runs[[paste0(s, "_start")]][runs$resident] <- NA_character_
    runs[[paste0(s, "_end")]][runs$resident] <- NA_character_
  }
}

# clean up
ebirdst_runs <- inner_join(runs, ebird_taxonomy, by = "species_code") %>%
  arrange(taxon_order) %>%
  mutate_at(vars(ends_with("start")), ~ ymd(paste0("2019-", .))) %>%
  select(run_name, species_code, scientific_name, common_name, resident,
         breeding_start, breeding_end,
         nonbreeding_start, nonbreeding_end,
         postbreeding_migration_start, postbreeding_migration_end,
         prebreeding_migration_start, prebreeding_migration_end,
         resident_start, resident_end)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
