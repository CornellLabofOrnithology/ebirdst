library(rnaturalearth)
library(raster)
library(sp)
devtools::load_all()

# Loading data
root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
species <- "yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83"
sp_path <- paste(root_path, species, sep = "")

# region for subsetting from rnaturalearth
wh_states <- ne_states(country = "United States of America")
wh_sub <- wh_states[wh_states$name == "Michigan", ]

mean_abund <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                   "_hr_2016_abundance_umean.tif"))
mean_occ <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                   "_hr_2016_occurrence_umean.tif"))
lower_abund <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                 "_hr_2016_abundance_lower.tif"))
upper_abund <- raster::stack(paste0(sp_path, "/results/tifs/", species,
                                 "_hr_2016_abundance_upper.tif"))

example_extent <- list(type = "polygon",
                  polygon = wh_sub)

example_mean_abund <- raster_st_subset(mean_abund, example_extent)
example_mean_occ <- raster_st_subset(mean_occ, example_extent)
example_lower_abund <- raster_st_subset(lower_abund, example_extent)
example_upper_abund <- raster_st_subset(upper_abund, example_extent)

writeRaster(example_mean_abund, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/tifs/", species, "_hr_2016_abundance_umean.tif"), format = "GTiff")

writeRaster(example_mean_occ, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/tifs/", species, "_hr_2016_occurrence_umean.tif"), format = "GTiff")

writeRaster(example_lower_abund, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/tifs/", species, "_hr_2016_abundance_lower.tif"), format = "GTiff")

writeRaster(example_upper_abund, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/tifs/", species, "_hr_2016_abundance_upper.tif"), format = "GTiff")

file.copy(paste0(sp_path, "/results/tifs/band_dates.csv"),
          paste0("/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/",
                 species, "/results/tifs/band_dates.csv"))

# abund_preds
# summary
summary <- read.csv(paste0(sp_path,
                           "/results/abund_preds/unpeeled_folds/summary.txt"),
                header = FALSE)
summary_spdf <- SpatialPointsDataFrame(summary[, c("V5", "V6")],
                                       summary,
                                       proj4string = CRS(proj4string(wh_sub)))
example_summary <- summary_spdf[!is.na(over(summary_spdf,
                                            as(wh_sub,
                                               "SpatialPolygons"))), ]@data
write.table(example_summary, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/abund_preds/summary.txt"), row.names = FALSE, col.names = FALSE,
  sep = ",", quote = FALSE)

# pis
pis <- load_pis(sp_path)

example_pis <- data_st_subset(loadpis, example_extent)
example_pis_names <- names(example_pis)[1:88]
example_pis$type <- "stixel"
example_pis$summary_type <- "pi"
example_pis <- example_pis[, c("type", "summary_type", example_pis_names)]

write.table(example_pis, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/abund_preds/pi.txt"), row.names = FALSE, col.names = FALSE,
  sep = ",", quote = FALSE)

# pds
pds <- load_pds(sp_path)

example_pds <- data_st_subset(pds, example_extent)

example_pds_names <- names(example_pds)[1:103]
example_pds$type <- "stixel"
example_pds$summary_type <- "pd"
example_pds <- example_pds[, c("type", "summary_type", example_pds_names)]

write.table(example_pds, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/abund_preds/pd.txt"), row.names = FALSE, col.names = FALSE,
  sep = ",", quote = FALSE)

# test
test <- read.csv(paste0(sp_path,
                        "/results/abund_preds/unpeeled_folds/test.pred.ave.txt"),
                 header = FALSE)
test_spdf <- SpatialPointsDataFrame(test[, c("V3", "V4")],
                                    test,
                                    proj4string = CRS(proj4string(wh_sub)))
example_test <- test_spdf[!is.na(over(test_spdf,
                                            as(wh_sub,
                                               "SpatialPolygons"))), ]@data
write.table(example_test, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/results/abund_preds/test.pred.ave.txt"), row.names = FALSE, col.names = FALSE,
  sep = ",", quote = FALSE)

# data pieces

file.copy(paste0(sp_path, "/data/", species, "_config.RData"),
          paste0("/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/",
                 species, "/data/", species, "_config.RData"))

file.copy(paste0(sp_path, "/data/", species, "_srd_raster_template.tif"),
          paste0("/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/",
                 species, "/data/", species, "_srd_raster_template.tif"))

rtest <- read.csv(paste0(sp_path, "/data/ebird.abund_", species,
                         "_erd.test.data.csv"),
                 header = TRUE)
rtest_spdf <- SpatialPointsDataFrame(rtest[, c("LONGITUDE", "LATITUDE")],
                                     rtest,
                                     proj4string = CRS(proj4string(wh_sub)))
example_rtest <- rtest_spdf[!is.na(over(rtest_spdf,
                                      as(wh_sub,
                                         "SpatialPolygons"))), ]@data
write.csv(example_rtest, paste0(
  "/Users/mta45/ebirdst/data-raw/clo-is-da-example-data/", species,
  "/data/ebird.abund_", species, "_erd.test.data.csv"), row.names = FALSE)
