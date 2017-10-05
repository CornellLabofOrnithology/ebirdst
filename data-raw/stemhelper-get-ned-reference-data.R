# gather the reference data
wh <- rnaturalearth::ne_countries(continent = c("North America",
                                                "South America"),
                                  scale = 50)
wh_states <- rnaturalearth::ne_states(iso_a2 = unique(wh@data$iso_a2))

mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

ned_wh_co_moll <- sp::spTransform(wh, mollweide)
ned_wh_st_moll <- sp::spTransform(wh_states, mollweide)

devtools::use_data(ned_wh_co_moll,
                   ned_wh_st_moll,
                   internal = TRUE)
