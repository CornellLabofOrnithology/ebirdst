library(xml2)
library(lubridate)
library(tidyverse)


runs_shorebird <- read_csv("data-raw/run_codes_shorebirds.csv",
                           col_names = "RUN_NAME") %>%
  mutate(SPECIES_CODE = str_extract(RUN_NAME, "[^-]+"))
runs_orig <- read_csv("data-raw/run_codes_2016.csv",
                      col_names = "RUN_NAME") %>%
  mutate(SPECIES_CODE = str_extract(RUN_NAME, "[^-]+"),
         SPECIES_CODE = if_else(SPECIES_CODE == "reevir1", "reevir",
                                SPECIES_CODE)) %>%
  filter(!SPECIES_CODE %in% runs_shorebird$SPECIES_CODE)
run_codes <- bind_rows(runs_shorebird, runs_orig) %>%
  arrange(SPECIES_CODE)

# attach taxonomic information
taxa <- read_csv("data-raw/taxonomy_2016.csv") %>%
  rename(scientific_name = SCI_NAME, common_name = PRIMARY_COM_NAME) %>%
  mutate_at(c("scientific_name", "common_name"),
            list(~ str_replace_all(., "_", " ")))
ebirdst_runs <- inner_join(run_codes, taxa, by = "SPECIES_CODE")
names(ebirdst_runs) <- names(ebirdst_runs) %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# add season definitions
# download regionSummaries-2018.csv from website manually
# https://ebird.org/science/status-and-trends/download-data/download?package=all-stats-regional
old_stats <- read_csv("data-raw/regionSummaries-2018.csv",
                      col_types = cols(.default = col_character())) %>%
  filter(!species_code %in% runs_shorebird$SPECIES_CODE) %>%
  mutate(species_code = if_else(species_code == "reevir1", "reevir",
                                species_code)) %>%
  select(species_code, season_name, start_dt, end_dt) %>%
  distinct() %>%
  mutate_at(c("start_dt", "end_dt"),
            funs(as.Date(str_replace(., "2015-", "2016-"))))
# shorebirds
shore_stats <- read_csv("data-raw/seasonDates-2018.csv",
                        col_types = cols(.default = col_character())) %>%
  mutate(start_dt = ymd(paste0("2016", start_week)),
         end_dt = ymd(paste0("2016", end_week))) %>%
  select(species_code, season_name, start_dt, end_dt) %>%
  distinct()
all_stats <- bind_rows(old_stats, shore_stats)
starts <- all_stats %>%
  select(species_code, season_name, start_dt) %>%
  mutate(season_name = paste0(season_name, "_start_dt")) %>%
  spread(key = season_name, value = start_dt)
ends <- all_stats %>%
  select(species_code, season_name, end_dt) %>%
  mutate(season_name = paste0(season_name, "_end_dt")) %>%
  spread(key = season_name, value = end_dt)
ebirdst_runs <- inner_join(starts, ends, by = "species_code") %>%
  select(species_code,
         breeding_start_dt, breeding_end_dt,
         nonbreeding_start_dt, nonbreeding_end_dt,
         postbreeding_migration_start_dt,
         postbreeding_migration_end_dt,
         prebreeding_migration_start_dt,
         prebreeding_migration_end_dt,
         year_round_start_dt, year_round_end_dt) %>%
  inner_join(ebirdst_runs, ., by = "species_code")

usethis::use_data(ebirdst_runs, overwrite = TRUE)
