library(tidyverse)
library(auk)
library(lubridate)
library(jsonlite)

pred_year <- file.path("data-raw", "config.json") %>%
  read_json(simplifyVector = TRUE) %>%
  pluck("SRD_PRED_YEAR")

runs <- read_csv("data-raw/seasons.csv") %>%
  rename_all(tolower) %>%
  filter(full_year_quality > 0) %>%
  rename(resident_quality = full_year_quality) %>%
  filter(!is.na(species_code))

# correctly na season dates
seasons <- c("breeding", "nonbreeding",
             "prebreeding_migration", "postbreeding_migration",
             "resident")
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

# default residents to full year
fy_resident <- runs$resident & runs$resident_quality > 0 &
  is.na(runs$resident_start) & is.na(runs$resident_end)
runs$resident_start[fy_resident] <- "01-04"
runs$resident_end[fy_resident] <- "12-28"

# clean up
ebirdst_runs <- runs %>%
  select(-common_name) %>%
  inner_join(ebird_taxonomy, by = "species_code") %>%
  arrange(taxon_order) %>%
  mutate_at(vars(ends_with("start")), ~ ymd(paste0(pred_year, "-", .))) %>%
  mutate_at(vars(ends_with("end")), ~ ymd(paste0(pred_year, "-", .))) %>%
  select(species_code, scientific_name, common_name, resident,
         breeding_quality, breeding_range_modeled,
         breeding_start, breeding_end,
         nonbreeding_quality, nonbreeding_range_modeled,
         nonbreeding_start, nonbreeding_end,
         postbreeding_migration_quality,
         postbreeding_migration_range_modeled,
         postbreeding_migration_start, postbreeding_migration_end,
         prebreeding_migration_quality,
         prebreeding_migration_range_modeled,
         prebreeding_migration_start, prebreeding_migration_end,
         resident_quality, resident_start, resident_end)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
