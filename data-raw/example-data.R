library(raster)
library(fs)
library(sf)
library(rnaturalearth)
library(tidyverse)
library(glue)
library(DBI)
library(RSQLite)
select <- dplyr::select

# source data
root_path <- path(rappdirs::user_data_dir("ebirdst"), "v2020")
species <- "yebsap"
sp_path <- path(root_path, species)

# destination
ex_species <- "yebsap_example"
ex_dir <- path(root_path, ex_species)
ex_cubes_dir <- path(ex_dir, "cubes")
ex_seasonal_dir <- path(ex_dir, "seasonal")
dir_create(ex_dir)
dir_create(ex_cubes_dir)
dir_create(ex_seasonal_dir)

# copy config file and databases
dir_ls(sp_path, regexp = "(json|db)$", recurse = FALSE, type = "file") %>%
  file_copy(ex_dir)

# raster template
r_tmplt <- path(sp_path, "srd_raster_template.tif") %>%
  raster()

# region for subsetting from rnaturalearth
state_mi_ll <- ne_states(iso_a2 = "US", returnclass = "sf") %>%
  filter(name == "Michigan") %>%
  st_geometry()
bb <- st_bbox(state_mi_ll)
state_mi <- state_mi_ll %>%
  st_transform(crs = projection(r_tmplt)) %>%
  as_Spatial()

# crop and mask rasters
tifs <- dir_ls(sp_path, regexp = "tif$", recurse = TRUE, type = "file")
for (f in tifs) {
  f_out <- str_replace_all(f, species, ex_species)
  r_ex <- f %>%
    stack() %>%
    crop(state_mi) %>%
    mask(state_mi) %>%
    writeRaster(filename = f_out, overwrite = TRUE,
                options = c("COMPRESS=DEFLATE","TILED=YES"))
}
dir_ls(path(sp_path, "cubes"), regexp = "csv$", recurse = FALSE, type = "file") %>%
  file_copy(ex_cubes_dir)

# database connections
pipd_db <- path(ex_dir, "stixel_summary.db")
pred_db <- path(ex_dir, "predictions.db")
pipd_con <- dbConnect(SQLite(), pipd_db)
pred_con <- dbConnect(SQLite(), pred_db)

# subset sqlite dbs to focal area
# predictions
sql <- str_glue("DELETE FROM predictions ",
                "WHERE longitude < {bb['xmin']} OR  longitude > {bb['xmax']} ",
                "OR latitude < {bb['ymin']} OR  latitude > {bb['ymax']};")
dbSendStatement(pred_con, sql)
# get sampling event ids that fall in polygon
sid <- tbl(pred_con, "predictions") %>%
  select(sampling_event_id, latitude, longitude) %>%
  collect() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_intersection(state_mi_ll) %>%
  st_drop_geometry()
copy_to(pred_con, sid)
# drop rows from other tables out of bounds of
sql <- str_glue("DELETE FROM predictions WHERE sampling_event_id NOT IN ",
                "(SELECT sampling_event_id FROM sid);")
dbSendStatement(pred_con, sql)
dbSendStatement(pred_con, "DROP TABLE sid;")
dbSendStatement(pred_con, "vacuum;")
dbDisconnect(pred_con)

# pipd
sql <- str_glue("DELETE FROM stixel_summary ",
                "WHERE lon < {bb['xmin']} OR  lon > {bb['xmax']} ",
                "OR lat < {bb['ymin']} OR  lat > {bb['ymax']};")
dbSendStatement(pipd_con, sql)
# get stixels centroids that fall in polygon
stixels <- tbl(pipd_con, "stixel_summary") %>%
  select(stixel_id, lat, lon) %>%
  collect() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_intersection(state_mi_ll) %>%
  st_drop_geometry()
copy_to(pipd_con, stixels)
# drop rows from other tables without out of bounds stixel centroids
for (t in c("occurrence_pds", "occurrence_pis",
            "count_pds", "count_pis",
            "stixel_summary")) {
  sql <- str_glue("DELETE FROM {t} WHERE stixel_id NOT IN ",
                  "(SELECT stixel_id FROM stixels);")
  dbSendStatement(pipd_con, sql)
}
dbSendStatement(pipd_con, "DROP TABLE stixels;")
dbSendStatement(pipd_con, "vacuum;")
dbDisconnect(pipd_con)


# copy to repo dir ---

repo_dir <- path("example-data", basename(ex_dir))
dir_copy(ex_dir, repo_dir)

# compress db files to meet GH size constrains
to_compress <- dir_ls(repo_dir, glob = "*.db")
wd <- getwd()
for (f in to_compress) {
  setwd(dirname(to_compress))
  zip(paste0(basename(f), ".zip"), basename(f))
  file_delete(basename(f))
  setwd(wd)
}

# file sizes must be below 100 mb
sizes <- dir_ls(repo_dir, type = "file", recurse = TRUE) %>%
  file_size() %>%
  as.numeric() %>%
  units::set_units("bytes")
stopifnot(sizes < units::set_units(90, "megabytes"))

# file list
file_list <- dir_ls(repo_dir, type = "file", recurse = TRUE) %>%
  str_remove("example-data/")
write_lines(file_list, "example-data/file-list.txt")
write_lines(file_list, "inst/extdata/example-data_file-list.txt")
