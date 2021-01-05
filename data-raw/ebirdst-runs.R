library(tidyverse)
library(auk)
library(lubridate)

runs <- read_csv("data-raw/seasons.csv") %>%
  rename_all(tolower) %>%
  filter(full_year_quality > 0) %>%
  rename(resident_quality = full_year_quality) %>%
  filter(!is.na(run_name))

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

# clean up
ebirdst_runs <- inner_join(runs, ebird_taxonomy, by = "species_code") %>%
  arrange(taxon_order) %>%
  mutate_at(vars(ends_with("start")), ~ ymd(paste0("2019-", .))) %>%
  select(run_name, species_code, scientific_name, common_name, resident,
         breeding_quality, breeding_start, breeding_end,
         nonbreeding_quality, nonbreeding_start, nonbreeding_end,
         postbreeding_migration_quality,
         postbreeding_migration_start, postbreeding_migration_end,
         prebreeding_migration_quality,
         prebreeding_migration_start, prebreeding_migration_end,
         resident_quality, resident_start, resident_end)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
