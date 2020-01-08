library(rnaturalearth)
library(sp)
library(sf)
library(data.table)
library(readr)
library(dplyr)

# source data
root_path <- rappdirs::user_data_dir("ebirdst")
species <- "yebsap-ERD2018-EBIRD_SCIENCE-20191030-3abe59ca"
sp_path <- file.path(root_path, species)

# destination
ex_species <- "yebsap-ERD2018-EBIRD_SCIENCE-20191030-3abe59ca-example"
ex_data_dir <- file.path(root_path, ex_species, "data")
ex_tif_dir <- file.path(root_path, ex_species, "results", "tifs")
ex_pred_dir <- file.path(root_path, ex_species, "results", "preds")
ex_stix_dir <- file.path(root_path, ex_species, "results", "stixels")
dir.create(ex_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ex_tif_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ex_pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ex_stix_dir, recursive = TRUE, showWarnings = FALSE)

# region for subsetting from rnaturalearth
us_states <- ne_states(country = "United States of America")
state_mi <- us_states[us_states$name == "Michigan", ]
state_mi_sf <- st_as_sf(state_mi) %>%
  st_transform(crs = 4326) %>%
  st_geometry()

# raster template
r_tmplt <- paste0(species, "_srd_raster_template.tif") %>%
  file.path(sp_path, "data", .) %>%
  raster::raster()
ex_tmplt <-  paste0(ex_species, "_srd_raster_template.tif") %>%
  file.path(ex_data_dir, .)
sp_ss <- spTransform(state_mi, raster::projection(r_tmplt))
r_tmplt_ex <- raster::crop(r_tmplt, sp_ss) %>%
  raster::mask(sp_ss) %>%
  raster::trim(values = NA) %>%
  raster::writeRaster(filename = ex_tmplt, overwrite = TRUE)

# subset tifs
f_tifs <- c("abundance_lower", "abundance_median", "abundance_upper",
            "count_median", "occurrence_median")
for (f in f_tifs) {
  r <- f %>%
    paste0(species, "_hr_2018_", ., ".tif") %>%
    file.path(sp_path, "results", "tifs", .) %>%
    raster::stack()
  f_ex <- f %>%
    paste0(ex_species, "_hr_2018_", ., ".tif") %>%
    file.path(ex_tif_dir, .)
  r_ex <- raster::crop(r, r_tmplt_ex) %>%
    raster::mask(r_tmplt_ex, filename = f_ex, overwrite = TRUE)
}
file.copy(file.path(sp_path, "results", "tifs", "band_dates.csv"),
          file.path(ex_tif_dir, "band_dates.csv"))

# subset function
subset_df <- function(f_in, f_out) {
  x <- fread(f_in)
  if ("LONGITUDE" %in% names(x)) {
    x_sf <- dplyr::select(x, LONGITUDE, LATITUDE)
    x_sf <- sf::st_as_sf(x, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  } else {
    x_sf <- dplyr::select(x, lon, lat)
    x_sf <- sf::st_as_sf(x, coords = c("lon", "lat"), crs = 4326)
  }
  is_in <- suppressMessages(
    sf::st_intersects(x_sf, state_mi_sf, sparse = FALSE)
  )
  x <- x[is_in[, 1, drop = TRUE], ]
  write_csv(x, f_out)
}
# test data
subset_df(file.path(sp_path, "data", paste0(species, "_test-data.csv")),
          file.path(ex_data_dir, paste0(ex_species, "_test-data.csv")))
# abundance predictions
subset_df(file.path(sp_path, "results", "preds", "test_pred_ave.txt"),
          file.path(ex_pred_dir, "test_pred_ave.txt"))
# summary
subset_df(file.path(sp_path, "results", "stixels", "summary.txt"),
          file.path(ex_stix_dir, "summary.txt"))
# pi
f_sum <- file.path(ex_stix_dir, "summary.txt") %>%
  fread() %>%
  select(stixel_id)
f_pi <- file.path(sp_path, "results", "stixels", "pi.txt") %>%
  fread()
inner_join(f_sum, f_pi, by = "stixel_id") %>%
  write_csv(file.path(ex_stix_dir, "pi.txt"))
