library(tidyverse)

run_names <- read_csv("data-raw/stem_erd2018_seasons.csv") %>%
  rename_all(tolower) %>%
  mutate(species_code = str_extract(run_name, "^.+(?=-ERD2018)"))
runs <- read_csv("data-raw/all-stats-regional-2019.csv")
# date fixes
fixes <- read_csv("data-raw/ebirdst-corrected-dates.csv")
for (i in seq.int(nrow(fixes))) {
  sp <- fixes$species_code[i]
  message(sp)
  season <- str_remove(fixes$season[i], "_end_dt|_start_dt")
  var <- str_extract(fixes$season[i], "end_dt|start_dt")
  dw <- fixes$date[i]
  dc <- fixes$correct_date[i]
  old_dt <- runs[runs$species_code == sp & runs$season_name == season, var]
  stopifnot(all(old_dt[[1]] == dw))
  runs[runs$species_code == sp & runs$season_name == season, var] <- dc
}
# correctly na season dates
seasons <- c("nonbreeding", "prebreeding_migration",
             "breeding", "postbreeding_migration")
for (i in seq.int(nrow(run_names))) {
  sp <- run_names$species_code[i]
  if (!sp %in% runs$species_code) {
    next()
  }
  for (s in seasons) {
    if (isTRUE(!run_names[[paste0(s, "_pass")]][i])) {
      message(paste(sp, s))
      data <- runs[runs$species_code == sp & runs$season_name == s,
                   "total_pop_percent", drop = TRUE]
      stopifnot(all(is.na(data)))
      runs[runs$species_code == sp & runs$season_name == s, "start_dt"] <- NA
      runs[runs$species_code == sp & runs$season_name == s, "end_dt"] <- NA
    }
  }
}
# write_csv(runs, "~/Desktop/all-stats-regional-2019.csv")

# continue
run_names <- select(run_names, species_code, run_name)
runs <- runs %>%
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
