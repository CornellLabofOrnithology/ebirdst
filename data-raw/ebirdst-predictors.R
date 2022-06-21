library(tidyverse)
library(jsonlite)

pl <- file.path("data-raw", "config.json") %>%
  read_json(simplifyVector = TRUE) %>%
  pluck("PREDICTOR_LIST") %>%
  ebirdst:::transform_predictor_names()

# weather labels
wv <- read_csv("data-raw/weather-variables.csv") %>%
  select(predictor = variable, label)

# predictors
ebirdst_predictors <- tibble(predictor = pl) %>%
  mutate(lc_class = str_remove_all(predictor, "_1500_[a-z]+$"),
         lc_class = if_else(str_detect(lc_class, "_fs_") |
                              str_detect(lc_class, "ntl") |
                              str_detect(lc_class, "gp_rtp"),
                            lc_class, NA_character_),
         lc_class = if_else(str_detect(lc_class, "ntl"),
                            "ntl", lc_class)) %>%
  # assign labels
  mutate(lc_class_label = case_when(
    # lccs landcover
    lc_class == "mcd12q1_lccs1_fs_c1" ~ "Barren",
    lc_class == "mcd12q1_lccs1_fs_c2" ~ "Permanent Snow and Ice",
    lc_class == "mcd12q1_lccs1_fs_c11" ~ "Evergreen Needleleaf Forests",
    lc_class == "mcd12q1_lccs1_fs_c12" ~ "Evergreen Broadleaf Forests",
    lc_class == "mcd12q1_lccs1_fs_c13" ~ "Deciduous Needleleaf Forests",
    lc_class == "mcd12q1_lccs1_fs_c14" ~ "Deciduous Broadleaf Forests",
    lc_class == "mcd12q1_lccs1_fs_c15" ~ "Mixed Broadleaf/Needleleaf Forests",
    lc_class == "mcd12q1_lccs1_fs_c16" ~ "Mixed Broadleaf Evergreen/Deciduous Forests",
    lc_class == "mcd12q1_lccs1_fs_c21" ~ "Open Forests",
    lc_class == "mcd12q1_lccs1_fs_c22" ~ "Sparse Forests",
    lc_class == "mcd12q1_lccs1_fs_c255" ~ "Unclassified",
    lc_class == "mcd12q1_lccs1_fs_c31" ~ "Dense Herbaceous",
    lc_class == "mcd12q1_lccs1_fs_c32" ~ "Sparse Herbaceous",
    lc_class == "mcd12q1_lccs1_fs_c41" ~ "Dense Shrublands",
    lc_class == "mcd12q1_lccs1_fs_c42" ~ "Shrubland/Grassland Mosaics",
    lc_class == "mcd12q1_lccs1_fs_c43" ~ "Sparse Shrublands",
    lc_class == "mcd12q1_lccs2_fs_c9" ~ "Urban and Built-up Lands",
    lc_class == "mcd12q1_lccs2_fs_c25" ~ "Forest/Cropland Mosaics",
    lc_class == "mcd12q1_lccs2_fs_c35" ~ "Natural Herbaceous/Croplands Mosaics",
    lc_class == "mcd12q1_lccs2_fs_c36" ~ "Herbaceous Croplands",
    lc_class == "mcd12q1_lccs3_fs_c27" ~ "Woody Wetlands",
    lc_class == "mcd12q1_lccs3_fs_c50" ~ "Herbaceous Wetlands",
    lc_class == "mcd12q1_lccs3_fs_c51" ~ "Tundra",
    # water cover
    lc_class == "astwbd_fs_c1" ~ "Ocean",
    lc_class == "astwbd_fs_c2" ~ "River",
    lc_class == "astwbd_fs_c3" ~ "Lakes",
    # roads
    lc_class == "gp_rtp_1" ~ "Highways (m/km2)",
    lc_class == "gp_rtp_2" ~ "Primary Roads (m/km2)",
    lc_class == "gp_rtp_3" ~ "Secondary Roads (m/km2)",
    lc_class == "gp_rtp_4" ~ "Tertiary Roads (m/km2)",
    lc_class == "gp_rtp_5" ~ "Local Roads (m/km2)",
    # intertidal
    lc_class == "intertidal_fs_c1" ~ "Tidal Mudflats",
    # ntl
    lc_class == "ntl" ~ "Nighttime Lights",
    TRUE ~ NA_character_)) %>%
  mutate(
    predictor_label = lc_class_label,
    predictor_label = coalesce(predictor_label,
                               str_replace_all(predictor, "_", " ") %>%
                                 str_to_title()),
    predictor_label = if_else(str_detect(predictor, "(ed|pland)$"),
                              paste(predictor_label,
                                    str_to_upper(str_extract(predictor, "(ed|pland)$"))),
                              predictor_label),
    predictor_label = str_replace(predictor_label, "Km", "(km)"),
    predictor_label = str_replace(predictor_label, "Hrs", "Hours"),
    predictor_label = str_replace(predictor_label, "Elev", "Elevation"),
    predictor_label = str_replace(predictor_label, "Sd", "SD"),
    predictor_label = ifelse(predictor == "cci", "CCI", predictor_label),
    predictor_label = ifelse(predictor == "solar_noon_diff_mid",
                             "Checklist Midpoint Time (Hours from Solar Noon)",
                             predictor_label),
    lc_class = if_else(str_detect(lc_class, "gp_rtp"),
                       "gp_rtp", lc_class),
    lc_class_label =  if_else(str_detect(lc_class, "gp_rtp"),
                              "Roads", lc_class_label),
    lc_class_label = dplyr::coalesce(lc_class_label, predictor_label)) %>%
  left_join(wv, by = "predictor") %>%
  mutate(predictor_label = coalesce(label, predictor_label),
         lc_class_label = coalesce(label, lc_class_label)) %>%
  select(predictor, predictor, predictor_label, lc_class, lc_class_label)

usethis::use_data(ebirdst_predictors, overwrite = TRUE)
