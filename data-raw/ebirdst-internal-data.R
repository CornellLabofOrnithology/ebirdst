library(dplyr)
library(stringr)
library(sf)
library(rnaturalearth)

# natural earth reference data
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
ne_scale <- 50

# countries with lakes removed, western hemisphere
# french guiana
fg <- ne_download(scale = ne_scale, category = "cultural",
                  type = "admin_0_map_units",
                  returnclass = "sf") %>%
  filter(ISO_A2 == "GF") %>%
  select(country_code = ISO_A2, country_name = NAME_EN)
# land border with lakes removed
ned_wh_co_moll <- ne_download(scale = ne_scale, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  #clean_names() %>%
  filter(CONTINENT %in% c("North America", "South America")) %>%
  select(country_code = ISO_A2, country_name = NAME_EN) %>%
  rbind(fg) %>%
  st_transform(crs = mollweide)

# states, north america
ned_wh_st_moll <- ne_download(scale = ne_scale, category = "cultural",
                        type = "admin_1_states_provinces_lakes",
                        returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN")) %>%
  select(country_code = iso_a2,
         state_code = iso_3166_2,
         state_name = gn_name) %>%
  st_transform(crs = mollweide)

# column names for tabular data
e <- new.env()
load("data-raw/yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83_config.RData",
     envir = e)

# summary
train_cov_means_names <- paste("train.cov.mean", e$PREDICTOR_LIST, sep = "_")
srd_cov_means_names <- paste("srd.cov.mean", e$PREDICTOR_LIST, sep = "_")
ebirdst_summary_names <- c("stixel", "data_type", "stixel.id", "srd.n",
                           "lon", "lat", "date",
                           "stixel_width", "stixel_height", "stixel_area",
                           "train.n",
                           "positive.ob_n",
                           "stixel_prevalence",
                           "mean_non_zero_count",
                           "binary_Kappa", "binary_AUC",
                           "binary.deviance_model", "binary.deviance_mean",
                           "binary.deviance_explained",
                           "pois.deviance_model", "pois.deviance_mean",
                           "posi.deviance_explained",
                           "total_EFFORT_HRS",
                           "total_EFFORT_DISTANCE_KM",
                           "total_NUMBER_OBSERVERS",
                           "train_elevation_mean",
                           train_cov_means_names, #k-covariate values
                           "train_covariate_entropy",
                           "srd_elevation_mean",
                           srd_cov_means_names, #k-covariate values
                           "srd_covariate_entropy",
                           "max_time")
ebirdst_summary_names <- ebirdst_summary_names %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# predictor importance
ebirdst_pi_names <- c("stixel", "data_type", "stixel.id", e$PI_VARS)
ebirdst_pi_names <- ebirdst_pi_names %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# partial dependence
ebirdst_pd_names <- c("stixel", "data_type", "stixel.id", "predictor",
                      paste0("y", 1:50), "spacer", paste0("x", 1:50))
ebirdst_pd_names <- ebirdst_pd_names %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# test data
ebirdst_td_names <- c("data_type", "row.id", "lon", "lat", "date", "obs",
                      "pi.mean", "pi.90", "pi.10", "pi.se", "pi.mu.mean",
                      "pi.mu.90", "pi.mu.10", "pi.mu.se", "pat", "pi.es")
ebirdst_td_names <- ebirdst_td_names %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# save as internal data files
usethis::use_data(ned_wh_co_moll, ned_wh_st_moll,
                  mollweide,
                  ebirdst_summary_names, ebirdst_pi_names,
                  ebirdst_pd_names, ebirdst_td_names,
                  internal = TRUE, overwrite = TRUE)
