library(reshape2)

run_codes <- read.csv("data-raw/run_codes_2016.csv", header = FALSE,
                      stringsAsFactors = FALSE)
names(run_codes) <- "RUN_NAME"

taxa <- read.csv("data-raw/taxonomy_2016.csv", stringsAsFactors = FALSE)

get_species_code <- function(x) {
  sc <- strsplit(x, "-")[[1]][1]

  if(sc == "reevir1") {
    sc <- "reevir"
  }

  return(sc)
}

run_codes$SPECIES_CODE <- apply(run_codes, 1, get_species_code)

runs_w_names <- merge(run_codes, taxa, by = "SPECIES_CODE")

runs_w_names$SCI_NAME <- apply(runs_w_names, 1, function(x) {
  gsub("_", " ", x["SCI_NAME"])
})

runs_w_names$PRIMARY_COM_NAME <- apply(runs_w_names, 1, function(x) {
  gsub("_", " ", x["PRIMARY_COM_NAME"])
})

# add seasons
all_stats <- read.csv("data-raw/regionSummaries-2018.csv", stringsAsFactors = FALSE)

stats_seasons <- all_stats[, c("species_code", "season_name", "start_dt",
                               "end_dt")]

stats_seasons$start_dt <- apply(stats_seasons, 1,
                                function(x){ gsub("2015-", "2016-",
                                                  x["start_dt"]) })
stats_seasons$end_dt <- apply(stats_seasons, 1,
                                function(x){ gsub("2015-", "2016-",
                                                  x["end_dt"]) })

ss_rs <- reshape(stats_seasons, idvar = "species_code", timevar = "season_name",
                 direction = "wide")

names(ss_rs) <- toupper(c("species_code", "breeding_start_dt", "breeding_end_dt",
                  "nonbreeding_start_dt", "nonbreeding_end_dt",
                  "postbreeding_migration_start_dt",
                  "postbreeding_migration_end_dt",
                  "prebreeding_migration_start_dt",
                  "prebreeding_migration_end_dt",
                  "year_round_start_dt", "year_round_end_dt"))

runs_w_names <- merge(runs_w_names, ss_rs, by = "SPECIES_CODE")

devtools::use_data(runs_w_names, overwrite = TRUE)
