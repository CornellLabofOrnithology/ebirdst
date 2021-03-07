library(dplyr)
library(stringr)
library(readr)
library(sf)
library(rnaturalearth)

# spatial ----

# natural earth reference data
prj_eck4 <- "+proj=eck4 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
prj_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
ne_scale <- 110

# project and recenter
recenter_sf <- function(x, crs) {
  stopifnot(inherits(x, c("sf", "sfc")))
  # find central meridian
  center <- as.numeric(stringr::str_extract(crs, "(?<=lon_0=)[-.0-9]+"))
  if (is.na(center) || length(center) != 1 || !is.numeric(center)) {
    stop("CRS has no central meridian term (+lon_0).")
  }

  # edge is offset from center by 180 degrees
  edge <- ifelse(center < 0, center + 180, center - 180)

  # create an very narrow sliver to clip out
  delta <- 1e-6
  clipper <- sf::st_bbox(c(xmin = edge - delta, xmax = edge + delta,
                           ymin = -90, ymax = 90),
                         crs = 4326)
  clipper <- sf::st_as_sfc(clipper)
  clipper <- suppressWarnings(smoothr::densify(clipper, max_distance = 1e-3))
  clipper <- sf::st_transform(clipper, crs = sf::st_crs(x))

  # cut then project
  x_proj <- sf::st_difference(x, clipper)
  sf::st_transform(x_proj, crs = crs)
}

# countries with lakes removed
ne_adm0_eck <- ne_download(scale = ne_scale, category = "cultural",
                            type = "admin_0_countries_lakes",
                            returnclass = "sf") %>%
  #clean_names() %>%
  select(country_code = ISO_A2, country_name = NAME_EN) %>%
  recenter_sf(crs = prj_eck4) %>%
  st_make_valid()

# states, north america
ne_adm1_eck <- ne_download(scale = ne_scale, category = "cultural",
                           type = "admin_1_states_provinces_lakes",
                           returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN")) %>%
  select(country_code = iso_a2,
         state_code = iso_3166_2,
         state_name = gn_name) %>%
  st_transform(crs = prj_eck4) %>%
  st_make_valid()


# colors ----

habitat_colors <- read_csv("data-raw/habitat-colors.csv")

# save ----
usethis::use_data(ne_adm0_eck, ne_adm1_eck, prj_eck4, prj_sinu,
                  habitat_colors,
                  internal = TRUE, overwrite = TRUE)
