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

# save as internal data files
usethis::use_data(ned_wh_co_moll, ned_wh_st_moll, mollweide,
                  internal = TRUE, overwrite = TRUE)
