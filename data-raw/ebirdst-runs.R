library(tidyverse)

run_names <- read_csv("data-raw/stem_erd2018_seasons.csv") %>%
  select(run_name = RUN_NAME) %>%
  mutate(species_code = str_extract(run_name, "^[^-]+"))
runs <- read_csv("data-raw/all-stats-regional-2019.csv") %>%
  inner_join(run_names, by = "species_code") %>%
  arrange(taxon_order) %>%
  distinct(run_name,
           species_code, scientific_name, common_name,
           season_name, start_dt, end_dt) %>%
  pivot_wider(names_from = season_name,
              values_from = c(start_dt, end_dt))
n <- names(runs)
n <- if_else(str_detect(n, "start_dt"),
             paste0(str_remove(n, "start_dt_"), "_start_dt"),
             n)
n <- if_else(str_detect(n, "end_dt"),
             paste0(str_remove(n, "end_dt_"), "_end_dt"),
             n)
names(runs) <- n
ebirdst_runs <- select(runs,
                       species_code, run_name, scientific_name, common_name,
                       breeding_start_dt, breeding_end_dt,
                       nonbreeding_start_dt, nonbreeding_end_dt,
                       postbreeding_migration_start_dt,
                       postbreeding_migration_end_dt,
                       prebreeding_migration_start_dt,
                       prebreeding_migration_end_dt,
                       year_round_start_dt, year_round_end_dt)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
