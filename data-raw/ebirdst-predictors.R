library(tidyverse)

e <- file.path("data-raw",
               "yebsap-ERD2018-EBIRD_SCIENCE-20191030-3abe59ca_config.rds") %>%
  readRDS()

# predictors
ebirdst_predictors <- tibble(predictor = e$PREDICTOR_LIST) %>%
  mutate(predictor_tidy = str_to_lower(predictor) %>%
           str_replace_all("\\.", "_"),
         lc_class = str_replace(predictor_tidy, "_1500_[a-z]+$", ""),
         lc_class = if_else(str_detect(lc_class, "_fs_") |
                              str_detect(lc_class, "ntl"),
                            lc_class, NA_character_),
         lc_class = if_else(str_detect(lc_class, "ntl"),
                            "ntl", lc_class)) %>%
  # assign labels
  dplyr::mutate(lc_class_label = dplyr::case_when(
    # umd landcover
    lc_class == "mcd12q1_umd_fs_c1" ~ "Evergreen Needleleaf Forests",
    lc_class == "mcd12q1_umd_fs_c10" ~ "Grasslands",
    lc_class == "mcd12q1_umd_fs_c11" ~ "Permanent Wetlands",
    lc_class == "mcd12q1_umd_fs_c12" ~ "Croplands",
    lc_class == "mcd12q1_umd_fs_c13" ~ "Urban and Built-up Lands",
    lc_class == "mcd12q1_umd_fs_c14" ~ "Cropland/Natural Vegetation Mosaics",
    lc_class == "mcd12q1_umd_fs_c15" ~ "Non-Vegetated Lands",
    lc_class == "mcd12q1_umd_fs_c2" ~ "Evergreen Broadleaf Forests",
    lc_class == "mcd12q1_umd_fs_c255" ~ "Unclassified",
    lc_class == "mcd12q1_umd_fs_c3" ~ "Deciduous Needleleaf Forests",
    lc_class == "mcd12q1_umd_fs_c4" ~ "Deciduous Broadleaf Forests",
    lc_class == "mcd12q1_umd_fs_c5" ~ "Mixed Forests",
    lc_class == "mcd12q1_umd_fs_c6" ~ "Closed Shrublands",
    lc_class == "mcd12q1_umd_fs_c7" ~ "Open Shrublands",
    lc_class == "mcd12q1_umd_fs_c8" ~ "Woody Savannas",
    lc_class == "mcd12q1_umd_fs_c9" ~ "Savannas",
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
    # esa landcovers
    lc_class == "esacci_lc_fs_c10" ~ "Cropland, rainfed",
    lc_class == "esacci_lc_fs_c100" ~ "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
    lc_class == "esacci_lc_fs_c11" ~ "Cropland, rainfed - Herbaceous cover",
    lc_class == "esacci_lc_fs_c110" ~ "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)",
    lc_class == "esacci_lc_fs_c12" ~ "Cropland, rainfed - Tree or shrub cover",
    lc_class == "esacci_lc_fs_c120" ~ "Shrubland",
    lc_class == "esacci_lc_fs_c121" ~ "Evergreen shrubland",
    lc_class == "esacci_lc_fs_c122" ~ "Deciduous shrubland",
    lc_class == "esacci_lc_fs_c130" ~ "Grassland",
    lc_class == "esacci_lc_fs_c140" ~ "Lichens and mosses",
    lc_class == "esacci_lc_fs_c150" ~ "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)",
    lc_class == "esacci_lc_fs_c152" ~ "Sparse shrub (<15%)",
    lc_class == "esacci_lc_fs_c153" ~ "Sparse herbaceous cover (<15%)",
    lc_class == "esacci_lc_fs_c160" ~ "Tree cover, flooded, fresh or brakish water",
    lc_class == "esacci_lc_fs_c170" ~ "Tree cover, flooded, saline water",
    lc_class == "esacci_lc_fs_c180" ~ "Shrub or herbaceous cover, flooded, fresh/saline/brakish water",
    lc_class == "esacci_lc_fs_c190" ~ "Urban areas",
    lc_class == "esacci_lc_fs_c20" ~ "Cropland, irrigated or post‐flooding",
    lc_class == "esacci_lc_fs_c200" ~ "Bare areas",
    lc_class == "esacci_lc_fs_c201" ~ "Consolidated bare areas",
    lc_class == "esacci_lc_fs_c202" ~ "Unconsolidated bare areas",
    lc_class == "esacci_lc_fs_c220" ~ "Permanent snow and ice",
    lc_class == "esacci_lc_fs_c30" ~ "Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)",
    lc_class == "esacci_lc_fs_c40" ~ "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)",
    lc_class == "esacci_lc_fs_c50" ~ "Tree cover, broadleaved, evergreen, closed to open (>15%)",
    lc_class == "esacci_lc_fs_c60" ~ "Tree cover, broadleaved, deciduous, closed to open (>15%)",
    lc_class == "esacci_lc_fs_c61" ~ "Tree cover, broadleaved, deciduous, closed (>40%)",
    lc_class == "esacci_lc_fs_c62" ~ "Tree cover, broadleaved, deciduous, open (15‐40%)",
    lc_class == "esacci_lc_fs_c70" ~ "Tree cover, needleleaved, evergreen, closed to open (>15%)",
    lc_class == "esacci_lc_fs_c71" ~ "Tree cover, needleleaved, evergreen, closed (>40%)",
    lc_class == "esacci_lc_fs_c72" ~ "Tree cover, needleleaved, evergreen, open (15‐40%)",
    lc_class == "esacci_lc_fs_c80" ~ "Tree cover, needleleaved, deciduous, closed to open (>15%)",
    lc_class == "esacci_lc_fs_c81" ~ "Tree cover, needleleaved, deciduous, closed (>40%)",
    lc_class == "esacci_lc_fs_c82" ~ "Tree cover, needleleaved, deciduous, open (15‐40%)",
    lc_class == "esacci_lc_fs_c90" ~ "Tree cover, mixed leaf type (broadleaved and needleleaved)",
    # water cover
    lc_class == "mod44w_oic_fs_c1" ~ "Ocean",
    lc_class == "mod44w_oic_fs_c2" ~ "Inland Water",
    lc_class == "mod44w_oic_fs_c3" ~ "Coastal Water",
    # intertidal
    lc_class == "intertidal_fs_c1" ~ "Tidal Mudflats",
    # ntl
    lc_class == "ntl" ~ "Nighttime Lights",
    TRUE ~ NA_character_)) %>%
  mutate(predictor_label = if_else(is.na(lc_class_label),
                                   str_replace_all(predictor_tidy, "_", " ") %>%
                                     str_to_title(),
                                   paste(lc_class_label,
                                         str_extract(predictor, "[A-Z]+$"))),
         predictor_label = str_replace(predictor_label, "Km", "(km)"),
         predictor_label = str_replace(predictor_label, "Hrs", "Hours"),
         predictor_label = str_replace(predictor_label, "Elev", "Elevation"),
         predictor_label = str_replace(predictor_label, "Sd", "SD"),
         predictor_label = str_replace(predictor_label, "Ntl", "Nighttime Lights"),
         lc_class_label = dplyr::coalesce(lc_class_label, predictor_label)) %>%
  select(predictor, predictor_tidy, predictor_label, lc_class, lc_class_label)

usethis::use_data(ebirdst_predictors, overwrite = TRUE)
