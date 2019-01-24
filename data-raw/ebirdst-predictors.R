library(dplyr)
library(stringr)

e <- new.env()
load("data-raw/yebsap-ERD2016-EBIRD_SCIENCE-20180729-7c8cec83_config.RData",
     envir = e)

# predictors
ebirdst_predictors <- tibble(predictor = e$PREDICTOR_LIST) %>%
  mutate(predictor_tidy = str_to_lower(predictor) %>%
           str_replace_all("\\.", "_"),
         lc_class = str_replace(predictor_tidy, "_1500_[a-z]+$", ""),
         lc_class = if_else(str_detect(lc_class, "_fs_"),
                            lc_class, NA_character_)) %>%
  # assign labels
  mutate(lc_class_label = case_when(
    # landcover
    lc_class == "umd_fs_c1" ~ "Evergreen Needleleaf Forest",
    lc_class == "umd_fs_c2" ~ "Evergreen Broadleaf Forest",
    lc_class == "umd_fs_c3" ~ "Deciduous Needleleaf Forest",
    lc_class == "umd_fs_c4" ~ "Deciduous Broadleaf Forest",
    lc_class == "umd_fs_c5" ~ "Mixed Forest",
    lc_class == "umd_fs_c6" ~ "Closed Shrublands",
    lc_class == "umd_fs_c7" ~ "Open Shrublands",
    lc_class == "umd_fs_c8" ~ "Woody Savannas",
    lc_class == "umd_fs_c9" ~ "Savannas",
    lc_class == "umd_fs_c10" ~ "Grasslands",
    lc_class == "umd_fs_c12" ~ "Croplands",
    lc_class == "umd_fs_c13" ~ "Urban",
    lc_class == "umd_fs_c16" ~ "Barren",
    # water cover
    lc_class == "modiswater_fs_c0" ~ "Shallow Ocean",
    lc_class == "modiswater_fs_c2" ~ "Ocean Coastlines and Lake Shores",
    lc_class == "modiswater_fs_c3" ~ "Shallow Inland Water",
    lc_class == "modiswater_fs_c5" ~ "Deep Inland Water",
    lc_class == "modiswater_fs_c6" ~ "Moderate Ocean",
    lc_class == "modiswater_fs_c7" ~ "Deep Ocean",
    TRUE ~ NA_character_)) %>%
  mutate(predictor_label = if_else(is.na(lc_class_label),
                                   str_replace_all(predictor_tidy, "_", " ") %>%
                                     str_to_title(),
                                   paste(lc_class_label,
                                         str_extract(predictor, "[A-Z]+$"))),
         predictor_label = str_replace(predictor_label, "Km", "(km)"),
         predictor_label = str_replace(predictor_label, "Hrs", "Hours"),
         predictor_label = str_replace(predictor_label, "Elev", "Elevation"),
         lc_class_label = coalesce(lc_class_label, predictor_label)) %>%
  select(predictor, predictor_tidy, predictor_label, lc_class, lc_class_label)
ebirdst_predictors <- as.data.frame(ebirdst_predictors,
                                    stringsAsFactors = FALSE)

usethis::use_data(ebirdst_predictors, overwrite = TRUE)
