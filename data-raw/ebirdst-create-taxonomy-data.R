library(tidyverse)

run_codes <- read_csv("data-raw/run_codes_2016.csv", col_names = "RUN_NAME") %>%
  mutate(SPECIES_CODE = str_extract(RUN_NAME, "[^-]+"),
         SPECIES_CODE = if_else(SPECIES_CODE == "reevir1", "reevir",
                                SPECIES_CODE))

# attach taxonomic information
taxa <- read_csv("data-raw/taxonomy_2016.csv") %>%
  rename(scientific_name = SCI_NAME, common_name = PRIMARY_COM_NAME) %>%
  mutate_at(c("scientific_name", "common_name"),
            funs(str_replace_all(., "_", " ")))
ebirdst_runs <- inner_join(run_codes, taxa, by = "SPECIES_CODE")
names(ebirdst_runs) <- names(ebirdst_runs) %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# add season definitions
all_stats <- read_csv("data-raw/regionSummaries-2018.csv",
                      col_types = cols(.default = col_character())) %>%
  select(species_code, season_name, start_dt, end_dt) %>%
  distinct() %>%
  mutate_at(c("start_dt", "end_dt"),
            funs(as.Date(str_replace(., "2015-", "2016-"))))
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
