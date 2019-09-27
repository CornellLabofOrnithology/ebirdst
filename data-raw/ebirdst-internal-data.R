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

# summary
ebirdst_summary_names <- c("stixel", "data_type", "stixel_id", "srd_n",
                           "lon", "lat", "date",
                           "stixel_width", "stixel_height", "stixel_area",
                           "train_n",
                           "positive_ob_n",
                           "stixel_prevalence",
                           "mean_non_zero_count",
                           "binary_kappa", "binary_auc",
                           "binary_deviance_model", "binary_deviance_mean",
                           "binary_deviance_explained",
                           "pois_deviance_model", "pois_deviance_mean",
                           "posi_deviance_explained",
                           "total_effort_hrs",
                           "total_effort_distance_km",
                           "total_number_observers",
                           "max_time")

# predictor importance
# pi column names
l <- readRDS("data-raw/woothr-ERD2018-BMEXP-20190726-4be38d37_config.rds")
ebirdst_pi_names <- c("stixel", "data_type", "response", "stixel_id",
                      l$PI_VARS, "encounter_rate")
ebirdst_pi_names <- ebirdst_pi_names %>%
  str_to_lower() %>%
  str_replace_all("\\.", "_")

# partial dependence
ebirdst_pd_names <- c("stixel", "data_type", "stixel_id", "predictor",
                      paste0("y", 1:50), "spacer", paste0("x", 1:50))

# test data
ebirdst_td_names <- c("data_type", "sampling_event_id",
                      "lon", "lat", "date", "obs",
                      "pi_mean", "pi_upper", "pi_lower", "pi_se",
                      "mu_mean", "mu_upper", "mu_lower", "mu_se",
                      "pi_mu_mean", "pi_mu_upper", "pi_mu_lower", "pi_mu_se",
                      "pat", "pi_es")

# save as internal data files
usethis::use_data(ned_wh_co_moll, ned_wh_st_moll,
                  mollweide,
                  ebirdst_summary_names, ebirdst_pi_names,
                  ebirdst_pd_names, ebirdst_td_names,
                  internal = TRUE, overwrite = TRUE)
