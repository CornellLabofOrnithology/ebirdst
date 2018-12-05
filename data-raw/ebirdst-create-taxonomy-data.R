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

devtools::use_data(runs_w_names, overwrite = TRUE)
