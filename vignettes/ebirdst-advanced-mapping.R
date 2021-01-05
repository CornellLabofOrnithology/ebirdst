## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "\\textwidth", 
                      fig.height = 5, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes local and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----libraries----------------------------------------------------------------
library(ebirdst)
library(raster)
library(sf)
library(smoothr)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# resolve namespace conflicts
select <- dplyr::select

## ----st-download--------------------------------------------------------------
# download to a temp directory for the vigette
# in practice, change to permanent directory the status and trends downloads
sp_path <- ebirdst_download(species = "example_data")
# load the abundance data
# this automaticaaly labels layers with their dates
abd <- load_raster("abundance", path = sp_path)

## ----ne-data, results="hide"--------------------------------------------------
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
ne_scale <- 50

# land polygon
ne_land <- ne_countries(scale = ne_scale, returnclass = "sf") %>%
  filter(continent %in% c("North America", "South America")) %>%
  st_set_precision(1e6) %>%
  st_union() %>% 
  st_geometry()
# function to subset other features to those  within this land area
wh_subset <- function(x) {
  in_wh <- as.logical(st_intersects(x, ne_land, sparse = FALSE))
  st_transform(x[in_wh], crs = mollweide)
}
# country lines
ne_country_lines <- ne_download(scale = ne_scale, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()
# state lines
ne_state_lines <- ne_download(scale = ne_scale, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()
# rivers
ne_rivers <- ne_download(scale = ne_scale, category = "physical",
                         type = "rivers_lake_centerlines",
                         returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()
# lakes
ne_lakes <- ne_download(scale = ne_scale, category = "physical",
                        type = "lakes",
                        returnclass = "sf") %>%
  st_geometry() %>%
  wh_subset()
ne_land <- st_transform(ne_land, crs = mollweide)

