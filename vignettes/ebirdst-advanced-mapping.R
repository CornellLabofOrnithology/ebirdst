## ----setup, include=FALSE------------------------------------------------
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

## ----libraries-----------------------------------------------------------
library(ebirdst)
library(raster)
library(velox)
library(sf)
library(smoothr)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# resolve namespace conflicts
select <- dplyr::select

## ----st-download---------------------------------------------------------
# download to a temp directory for the vigette
# in practice, change to permanent directory the status and trends downloads
sp_path <- download_data(species = "example_data")
# load the abundance data
# this automaticaaly labels layers with their dates
abd <- load_raster("abundance_umean", path = sp_path)

## ----ne-data, results="hide"---------------------------------------------
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

## ----season-defs---------------------------------------------------------
# subset to the yellow-bellied sapsucker season definitions
yebsap_dates <- filter(ebirdst_runs, species_code == "yebsap") %>% 
  # just keep the seasonal definition columns
  select(setdiff(matches("(start)|(end)"), matches("year_round"))) %>% 
  # transpose
  gather("label", "date") %>% 
  # spread data so start and end dates are in separate columns
  separate(label, c("season", "start_end"), "_(?=s|e)") %>% 
  spread(start_end, date) %>% 
  select(season, start_dt, end_dt)
# did the season pass review
yebsap_dates <- mutate(yebsap_dates, pass = !(is.na(start_dt) | is.na(end_dt)))
yebsap_dates

## ----season-assignment---------------------------------------------------
# dates for each abundance layer
weeks <- parse_raster_dates(abd)
# assign to seasons
weeks_season <- rep(NA_character_, length(weeks))
for (i in seq_len(nrow(yebsap_dates))) {
  s <- yebsap_dates[i, ]
  # skip seasona assignment if season failed
  if (!s$pass) {
    next()
  }
  # handle seasons cross jan 1 separately
  if (s$start_dt <= s$end_dt) {
    in_season <- weeks >= s$start_dt & weeks <= s$end_dt
  } else {
    in_season <- weeks >= s$start_dt | weeks <= s$end_dt
  }
  weeks_season[in_season] <- s$season
}
table(weeks_season)

## ----seasonal-average----------------------------------------------------
# drop weeks not assigned to season
week_pass <- !is.na(weeks_season)
abd <- abd[[which(week_pass)]]
weeks <- weeks[week_pass]
weeks_season <- weeks_season[week_pass]
# average over weeks in season
mean_season <- function(s) {
  calc(abd[[which(weeks_season == s)]], mean, na.rm = TRUE)
}
seasons <- unique(weeks_season)
abd_season <- lapply(seasons, mean_season) %>% 
  stack() %>% 
  setNames(seasons)
abd_season

## ----split-seasons-------------------------------------------------------
migration_threshold <- 0.4
mig_seasons <- c("prebreeding_migration", "postbreeding_migration")
if (all(mig_seasons %in% names(abd_season))) {
  # identify areas with abundance in only one season
  abd_nz <- abd_season[[mig_seasons]] > 0
  just_pre <- mask(abd_nz[["prebreeding_migration"]],
                   abd_nz[["postbreeding_migration"]], 
                   maskvalue = 1)
  just_post <- mask(abd_nz[["postbreeding_migration"]],
                   abd_nz[["prebreeding_migration"]], 
                   maskvalue = 1)
  # count the number of cells with abundance in only one season
  n_just <- cellStats(stack(just_pre, just_post), sum)
  n_all <- cellStats(abd_nz, sum)
  # is the proportion of one season cells above the 40% threshold
  split_migration <- max(n_just / n_all, na.rm = TRUE) >= migration_threshold
} else {
  split_migration <- FALSE
}
n_just / n_all
split_migration

## ----show-yr-------------------------------------------------------------
threshold_yearround <- 0.01
# decide whether to show year-round layer
if (nlayers(abd_season) == 4) {
  # annual abundance
  abd_yr <- calc(abd, fun = mean, na.rm = TRUE)
  # mask out cells that aren't occupied year-round
  year_round <- calc(abd_season > 0, fun = sum, na.rm = TRUE) == 4
  abd_yr_mask <- mask(abd_yr, year_round, maskvalue = 0)
  # determine proportion of celss that are occupied year round
  n_yr <- cellStats(abd_yr_mask > 0, sum)
  n_an <- cellStats(abd_yr > 0, sum)
  # only show year round abundance if it's above 1% of range threshold
  show_yearround <- ((n_yr / n_an) >= threshold_yearround)
} else {
  show_yearround <- FALSE
}
show_yearround

## ----breaks--------------------------------------------------------------
bin_breaks <- calc_bins(abd_season)

## ----abd-map-------------------------------------------------------------
# project the abundance data to mollweide
# use nearest neighbour resampling to preserve true zeros
abd_season_proj <- projectRaster(abd_season, crs = mollweide, method = "ngb")
# determine spatial extent for plotting
ext <- calc_full_extent(abd_season_proj)
# set the plotting order of the seasons
season_order <- c("postbreeding_migration", "prebreeding_migration", 
                  "nonbreeding", "breeding")

# prediction region, cells with predicted value in at least one week
pred_region <- calc(abd_season_proj, mean, na.rm = TRUE)
# mask to land area
ne_land_buffer <- st_buffer(ne_land, dist = max(res(pred_region)) / 2)
pred_region <- mask(pred_region, as_Spatial(ne_land_buffer))

# remove zeros from abundnace layers
abd_no_zero <- subs(abd_season_proj, data.frame(from = 0, to = NA), 
                    subsWithNA = FALSE)

# set up plot area
par(mar = c(0 , 0, 0, 0))
plot(ne_land, col = "#eeeeee", border = NA, 
     xlim = c(ext@xmin, ext@xmax),
     ylim = c(ext@ymin, ext@ymax))
# prediction region and explicit zeros
plot(pred_region, col = "#dddddd", maxpixels = raster::ncell(pred_region),
     legend = FALSE, add = TRUE)
# lakes
plot(ne_lakes, col = "#ffffff", border =  "#444444", lwd = 0.5, add = TRUE)
# land border
plot(ne_land, col = NA, border = "#444444", lwd = 0.5, add = TRUE)
# seasonal layer
plot_seasons <- intersect(season_order, names(abd_no_zero))
for (s in plot_seasons) {
  # handle splitting of migration seasons into different colors
  if (!split_migration && s %in% c("prebreeding_migration", 
                                   "postbreeding_migration")) {
    pal_season <- "migration"
    
  } else {
    pal_season <- s
  }
  pal <- abundance_palette(length(bin_breaks$bins) - 1, pal_season)
  plot(abd_no_zero[[s]], col = pal, breaks = bin_breaks$bins,
       maxpixels = ncell(abd_no_zero[[s]]),
       legend = FALSE, add = TRUE)
}
# year round
if (show_yearround) {
  year_round_proj <- projectRaster(year_round, crs = mollweide, method = "ngb")
  plot(year_round_proj, 
       col = abundance_palette(length(bin_breaks$bins) - 1, "year_round"), 
       breaks = bin_breaks$bins,
       maxpixels = ncell(year_round_proj),
       legend = FALSE, add = TRUE)
}
# linework
plot(ne_rivers, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 2, add = TRUE)

# legends
legend_seasons <- plot_seasons
if (split_migration) {
  legend_seasons[legend_seasons %in% c("prebreeding_migration", 
                                       "postbreeding_migration")] <- "migration"
  legend_seasons <- unique(legend_seasons)
}
if (show_yearround) {
  legend_seasons <- c(legend_seasons, "year_round")
}
# thin out labels
lbl_at <- bin_breaks$bins^bin_breaks$power
lbl_at <- c(min(lbl_at), median(lbl_at), max(lbl_at))
lbl <- lbl_at^(1 / bin_breaks$power)
lbl <- format(round(lbl, 2), nsmall = 2)
# plot legends
for (i in seq_along(legend_seasons)) {
  pal <- abundance_palette(length(bin_breaks$bins) - 1, legend_seasons[i])
  if (i == 1) {
    axis_args <- list(at = lbl_at, labels = lbl, line = -1,
                      cex.axis = 0.75, lwd = 0)
  } else {
    axis_args <- list(at = lbl_at, labels = rep("", 3),
                      cex.axis = 0.75, lwd = 0)
  }
  legend_title <- legend_seasons[i] %>% 
    str_replace_all("_", " ") %>% 
    str_to_title()
  fields::image.plot(zlim = range(bin_breaks$bins^bin_breaks$power), 
                   legend.only = TRUE, 
                   breaks = bin_breaks$bins^bin_breaks$power, col = pal,
                   smallplot = c(0.05, 0.35, 0.01 + 0.06 * i, 0.03 + 0.06 * i),
                   horizontal = TRUE,
                   axis.args = axis_args,
                   legend.args = list(text = legend_title, side = 3, 
                                      cex = 0.9, col = "black", line = 0.1))
}
title("Yellow-bellied Sapsucker Relative Abundance (birds per km/hr)", 
      line = -1, cex.main = 1)

## ----raster-to-polygon---------------------------------------------------
# aggregate
abd_season_agg <- aggregate(abd_season_proj, fact = 3)
# raster to polygon, one season at a time
range <- list()
pred_area <- list()
for (s in names(abd_season_agg)) {
  # range
  range[[s]] <- rasterToPolygons(abd_season_agg[[s]], 
                                 fun = function(y) {y > 0}, 
                                 digits = 6) %>% 
    st_as_sfc() %>% 
    # combine polygon pieces into a single multipolygon
    st_set_precision(1e6) %>% 
    st_union() %>% 
    st_sf() %>% 
    # tag layers with season
    mutate(season = s, layer = "range")
  # prediction area
  pred_area[[s]] <- rasterToPolygons(abd_season_agg[[s]], 
                                     fun = function(y) {!is.na(y)}, 
                                     digits = 6) %>% 
    st_as_sfc() %>% 
    # combine polygon pieces into a single multipolygon
    st_set_precision(1e6) %>% 
    st_union() %>% 
    st_sf() %>% 
    # tag layers with season
    mutate(season = s, layer = "prediction_area")
}
# combine the sf objects for all seasons
range <- rbind(do.call(rbind, range), do.call(rbind, pred_area))
row.names(range) <- NULL
print(range)

## ----smoothr-------------------------------------------------------------
# clean and smooth
cell_area <- (1.5 * prod(res(abd_season_agg)))
range_smooth <- range %>% 
  # drop fragment polygons smaller than 1.5 times the aggregated cell size
  drop_crumbs(threshold = cell_area) %>% 
  # drop holes in polygons smaller than 1.5 times the aggregated cell size
  fill_holes(threshold = cell_area) %>% 
  # smooth the polygon edges
  smooth(method = "ksmooth", smoothness = 2)
# clip zeros to land border, range to buffered land to handle coastal species
range_split <- split(range_smooth, range_smooth$layer)
range_smooth <- rbind(
  st_intersection(range_split$range, ne_land_buffer),
  st_intersection(range_split$prediction_area, ne_land))

## ----range-map-----------------------------------------------------------
# range map color palette
range_palette <- c(nonbreeding = "#1d6996",
                   prebreeding_migration = "#73af48",
                   breeding = "#cc503e",
                   postbreeding_migration = "#edad08",
                   migration = "#edad08",
                   year_round = "#6f4070")

# set up plot area
par(mar = c(0 , 0, 0, 0))
plot(ne_land, col = "#eeeeee", border = NA, 
     xlim = c(ext@xmin, ext@xmax),
     ylim = c(ext@ymin, ext@ymax))
# prediction region and explicit zeros
annual_pred_area <- filter(range_smooth, layer == "prediction_area") %>% 
  st_union()
plot(annual_pred_area, col = "#dddddd", border = NA, add = TRUE)
# lakes
plot(ne_lakes, col = "#ffffff", border =  "#444444", lwd = 0.5, add = TRUE)
# land border
plot(ne_land, col = NA, border = "#444444", lwd = 0.5, add = TRUE)
# seasonal layer
for (s in intersect(season_order, unique(range_smooth$season))) {
  # handle splitting of migration seasons into different colors
  if (!split_migration && s %in% c("prebreeding_migration", 
                                   "postbreeding_migration")) {
    col_season <- "migration"
  } else {
    col_season <- s
  }
  rng_season <- filter(range_smooth, season == s, layer == "range") %>% 
    st_geometry()
  plot(rng_season, col = range_palette[col_season], border = NA, add = TRUE)
}
# year round
if (show_yearround) {
  # find common area between all seasons
  range_combined <- filter(range_smooth, layer == "range")
  range_yearround <- range_combined[1, ]
  range_combined <- sf::st_geometry(range_combined)
  for (i in 2:length(range_combined)) {
    range_yearround <- sf::st_intersection(range_yearround, range_combined[i])
  }
  plot(st_geometry(range_yearround), 
       col = range_palette["year_round"], border = NA, 
       add = TRUE)
}
# linework
plot(ne_rivers, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 2, add = TRUE)

# legend
rng_legend <- rev(range_palette[legend_seasons])
names(rng_legend) <- names(rng_legend) %>% 
    str_replace_all("_", " ") %>% 
    str_to_title()
legend("bottomleft", legend = names(rng_legend), fill = rng_legend)
title("Yellow-bellied Sapsucker Seasonal Range Map", 
      line = -1, cex.main = 1)

## ----counties------------------------------------------------------------
mi_counties <- getData("GADM", country = "USA", level = 2, path = tempdir()) %>% 
  st_as_sf() %>% 
  filter(NAME_1 == "Michigan") %>% 
  select(county = NAME_2, county_code = HASC_2) %>% 
  st_transform(crs = projection(abd)) %>% 
  # remove lakes which aren't true counties
  filter(county_code != "US.MI.WB")

## ----aggregate-----------------------------------------------------------
abd_season_agg <- aggregate(abd_season, fact = 3, fun = mean, na.rm = TRUE)
abd_agg <- aggregate(abd, fact = 3, fun = mean, na.rm = TRUE)

## ----mean-rel-abd--------------------------------------------------------
vx <- velox(abd_season_agg)
abd_mean_region <- vx$extract(mi_counties, 
                              fun = function(x) {mean(x, na.rm = TRUE)}, 
                              df = TRUE) %>% 
  # assign season names
  setNames(c("id", names(abd_season_agg))) %>% 
  # attach county attributes
  mutate(county_code = mi_counties$county_code) %>% 
  select(-id)
# transpose to one season per row
abd_mean_region <- abd_mean_region %>% 
  gather(season, mean_abundance, -county_code)
head(abd_mean_region)

## ----mean-rel-abd-map----------------------------------------------------
# join back to county boundaries
abd_mean_proj <- abd_mean_region %>% 
  inner_join(mi_counties, ., by = "county_code") %>% 
  # transform to mollweide for plotting
  st_transform(crs = mollweide)

# plot
ggplot(abd_mean_proj %>% filter(season %in% c("breeding", "nonbreeding"))) +
  # abundance data
  geom_sf(aes(fill = mean_abundance)) +
  scale_fill_viridis_c(trans = "sqrt") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
  facet_wrap(~ season, ncol = 2) +
  labs(title = "Seasonal mean relative abundance in MI counties",
       fill = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "bottom")

## ----week-occ------------------------------------------------------------
threshold_occupied <- 0.05

# weekly percent occupied
calc_week_occupied <- function(wk) {
  vx <- velox(abd_agg[[wk]])
  vx$extract(mi_counties, 
             fun = function(x) {sum(x > 0, na.rm = TRUE) / length(x)}, 
             df = TRUE) %>%
    # assign week, season, and county
    transmute(week_date = weeks[wk],
              season = weeks_season[wk],
              county_code = mi_counties$county_code,
              pct_occupied = coalesce(out, 0))
}
# process weeks independently otherwise we'll run into memory issues
week_occupied <- lapply(seq_len(nlayers(abd_agg)), calc_week_occupied) %>% 
  bind_rows() %>% 
  mutate(occupied = (pct_occupied > threshold_occupied))

## ----day-occ-------------------------------------------------------------
days_occupation <- week_occupied %>%
  group_by(county_code, season) %>%
  summarise(weeks_occ = sum(occupied)) %>%
  ungroup() %>%
  mutate(days_occ = round(7 * weeks_occ)) %>% 
  select(county_code, season, days_occ)
head(days_occupation)

## ----day-occ-map---------------------------------------------------------
# join back to county boundaries
days_occ_proj <- days_occupation %>% 
  inner_join(mi_counties, ., by = "county_code") %>% 
  # transform to mollweide for plotting
  st_transform(crs = mollweide)

# plot
ggplot(days_occ_proj %>% filter(season %in% c("breeding", "nonbreeding"))) +
  # abundance data
  geom_sf(aes(fill = days_occ)) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
  facet_wrap(~ season, ncol = 2) +
  labs(title = "Seasonal days of occupation in MI counties",
       fill = "Days of occupation") +
  theme_bw() +
  theme(legend.position = "bottom")

