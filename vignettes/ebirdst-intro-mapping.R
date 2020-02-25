## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "\\textwidth", 
                      fig.height = 4, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes locally and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----load_raster_stack--------------------------------------------------------
library(ebirdst)
library(raster)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(viridisLite)
# handle namespace conflicts
extract <- raster::extract

# Currently, example data is available on a public s3 bucket. The following 
# ebirdst_download() function copies the species results to a selected path and 
# returns the full path of the results. In this vignette, we'll use example data 
# for Yellow-bellied Sapsucker
sp_path <- ebirdst_download(species = "example_data")

# load trimmed mean abundances and label dates
abunds <- load_raster("abundance", path = sp_path)

# use parse_raster_dates() to get actual date objects for each layer
date_vector <- parse_raster_dates(abunds)
print(date_vector)

## ----aggregate----------------------------------------------------------------
abunds_low_res <- aggregate(abunds, fact = 3, fun = mean)
# native resolution
res(abunds)
# resolution after aggregation
res(abunds_low_res)

## ----project_stack------------------------------------------------------------
# define mollweide projection
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

# project single layer from stack to mollweide
week38_moll <- projectRaster(abunds[[38]], crs = mollweide, method = "ngb")
week38_moll <- trim(week38_moll)

# optionally, you can project an entire stack, but it takes much longer
# abund_moll <- projectRaster(abund, crs = mollweide, method = "ngb")

# map single layer with full annual extent
par(mar = c(0, 0, 0, 2))
plot(week38_moll, 
     col = abundance_palette(10, season = "weekly"), 
     axes = FALSE, box = FALSE, 
     maxpixels = ncell(week38_moll))

## ----map_occurrence-----------------------------------------------------------
occs <- load_raster("occurrence", path = sp_path)

# select a week in the summer
occ <- occs[[26]]

# create breaks every 0.1 from 0 to 1
occ_bins <- seq(0, 1, by = 0.1)
occ_moll <- projectRaster(occ, crs = mollweide, method = "ngb")
occ_moll <- trim(occ_moll)

par(mar = c(0, 0, 0, 2), cex = 0.9)
plot(occ_moll, 
     breaks = occ_bins, 
     col = abundance_palette(length(occ_bins) - 1, season = "weekly"),
     axes = FALSE, box = FALSE, 
     maxpixels = ncell(occ_moll),
     legend.width = 2, legend.shrink = 0.97)

## ----map_linear---------------------------------------------------------------
year_max <- max(maxValue(abunds), na.rm = TRUE)

week14_moll <- projectRaster(abunds[[14]], crs = mollweide, method = "ngb")
week14_moll <- trim(week14_moll)

# set graphical params
par(mfrow = c(1, 2), mar = c(0, 0, 0, 4))

# use raster bounding box to set the spatial extent for the plot
bb <- st_as_sfc(st_bbox(trim(week14_moll)))
plot(bb, col = "white", border = "white")
# plot the abundance
plot(week38_moll, zlim = c(0, year_max), 
     col = abundance_palette(20, season = "weekly"), 
     maxpixels = ncell(week38_moll),
     axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)

# do the same for week 14
par(mar = c(0, 0, 0, 4))
bb <- st_as_sfc(st_bbox(trim(week14_moll)))
plot(bb, col = "white", border = "white")
plot(week14_moll, zlim = c(0, year_max),
     col = abundance_palette(20, season = "weekly"), 
     maxpixels = ncell(week14_moll),
     axes = FALSE, box = FALSE, legend.shrink = 0.97, add = TRUE)

## ----map_bins-----------------------------------------------------------------
# calculate ideal color bins for abundance values
year_bins <- calc_bins(abunds)

# plot
par(mfrow = c(1, 2), mar = c(0, 0, 0, 6))
plot(st_as_sfc(st_bbox(week38_moll)), col = "white", border = "white")
plot(week38_moll, 
     breaks = year_bins$bins, 
     col = abundance_palette(length(year_bins$bins) - 1, season = "weekly"), 
     maxpixels = ncell(week38_moll), 
     axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)
par(mar = c(0, 0, 0, 6))
plot(st_as_sfc(st_bbox(trim(week14_moll))), col = "white", border = "white")
plot(week14_moll, 
     breaks = year_bins$bins, 
     col = abundance_palette(length(year_bins$bins) - 1, season = "weekly"), 
     maxpixels = ncell(week14_moll), 
     axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)

# create a thinner set of labels
bin_labels <- format(round(year_bins$bins, 2), nsmall = 2)
bin_labels[!(bin_labels %in% c(bin_labels[1],
                               bin_labels[round((length(bin_labels) / 2)) + 1],
                               bin_labels[length(bin_labels)]))] <- ""

# plot legend
plot(week38_moll^year_bins$power, legend.only = TRUE, 
     col = abundance_palette(length(year_bins$bins) - 1, season = "weekly"), 
     breaks = year_bins$bins ^ year_bins$power, lab.breaks = bin_labels,
     legend.shrink = 0.97, legend.width = 2,  
     axis.args = list(cex.axis = 0.9, lwd.ticks = 0))

## ----map_w_es, out.width = NULL-----------------------------------------------
# to add context, let's pull in some reference data to add
wh_states <- ne_states(country = c("United States of America", "Canada"),
                    returnclass = "sf") %>% 
  st_transform(crs = mollweide) %>% 
  st_geometry()

# we'll plot a week in the middle of summer
week26_moll <- projectRaster(abunds[[26]], crs = mollweide, method = "ngb")
week26_moll <- trim(week26_moll)

# set graphics params
par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))

# use the extent object to set the spatial extent for the plot
plot(st_as_sfc(st_bbox(week26_moll)), col = "white", border = "white")

# add background spatial context
plot(wh_states, col = "#eeeeee", border = NA, add = TRUE)

# plot zeroes as gray
plot(week26_moll == 0, col = "#dddddd", 
     maxpixels = ncell(week26_moll),
     axes = FALSE, legend = FALSE, add = TRUE)

# define color bins
qcol <- abundance_palette(length(year_bins$bins) - 1, "weekly")

# plot abundances
plot(week26_moll, breaks = year_bins$bins, col = qcol, 
     maxpixels = ncell(week26_moll),
     axes = FALSE, legend = FALSE, add = TRUE)

# for legend, create a smaller set of bin labels
bin_labels <- format(round(year_bins$bins, 2), nsmall = 2)
bin_labels[!(bin_labels %in% c(bin_labels[1],
                               bin_labels[round((length(bin_labels) / 2)) + 1],
                               bin_labels[length(bin_labels)]))] <- ""
bin_labels <- c("0", bin_labels)

# create colors that include gray for 0
lcol <- c("#dddddd", qcol)

# set legend such that color ramp appears linearly
ltq <- seq(from = year_bins$bins[1], 
           to = year_bins$bins[length(year_bins$bins)],
           length.out = length(year_bins$bins))
ltq <- c(0, ltq)

# plot legend
plot(week26_moll^year_bins$power, legend.only = TRUE,
     col = lcol, breaks = ltq^year_bins$power, 
     lab.breaks = bin_labels, 
     legend.shrink = 0.97, legend.width = 2, 
     axis.args = list(cex.axis = 0.9, lwd.ticks = 0))

# add state boundaries on top
plot(wh_states, border = "white", lwd = 1.5, add = TRUE)

## ----map_confidence_band------------------------------------------------------
# load lower and upper stacks
# load trimmed mean abundances and label dates
lower <- load_raster("abundance_lower", path = sp_path)
upper <- load_raster("abundance_upper", path = sp_path)

# calculate band width
conf_band <- upper[[26]] - lower[[26]]

conf_week26 <- projectRaster(conf_band, crs = mollweide, method = "ngb")

par(mar = c(0, 0, 0, 2))
plot(trim(conf_week26), col = magma(20), 
     maxpixel = ncell(conf_week26),
     axes = FALSE, box = FALSE)

## ----trajectories-------------------------------------------------------------
# set a point
pt <- st_point(c(-88.1, 46.7)) %>% 
  st_sfc(crs = 4326) %>% 
  st_transform(crs = projection(abunds)) %>% 
  st_coordinates()

# extract
abund_traj <- extract(abunds, pt, fun = mean, na.rm = TRUE)[1, ]
upper_traj <- extract(upper, pt, fun = mean, na.rm = TRUE)[1, ]
lower_traj <- extract(lower, pt, fun = mean, na.rm = TRUE)[1, ]

# plot trajectories
plot_frame <- data.frame(x = 1:length(abund_traj),
                         y = unname(abund_traj),
                         upper = unname(upper_traj),
                         lower = unname(lower_traj))

ggplot(plot_frame, aes(x, y)) +
  geom_line(data = plot_frame) +
  geom_ribbon(data = plot_frame, 
              aes(ymin = lower, ymax = upper), 
              alpha = 0.3) +
  ylab("Expected Count") +
  xlab("Week") +
  theme_light()

## ----trajectories_region------------------------------------------------------
# set an extent based on polygon
mi <- ne_states(country = "United States of America", returnclass = "sf") %>% 
  st_transform(crs = projection(abunds))
mi <- mi[mi$name == "Michigan"]

# extract
# because we're using a region, we get lots of values that we need to average
abund_traj <- extract(abunds, mi, fun = mean, na.rm = TRUE)
abund_traj <- apply(abund_traj, 2, mean, na.rm = TRUE)

upper_traj <- extract(upper, mi, fun = mean, na.rm = TRUE)
upper_traj <- apply(upper_traj, 2, mean, na.rm = TRUE)

lower_traj <- extract(lower, mi, fun = mean, na.rm = TRUE)
lower_traj <- apply(lower_traj, 2, mean, na.rm = TRUE)

# plot trajectories
plot_frame <- data.frame(x = 1:length(abund_traj),
                         y = unname(abund_traj),
                         upper = unname(upper_traj),
                         lower = unname(lower_traj))

ggplot(plot_frame, aes(x, y)) +
  geom_line(data = plot_frame) +
  geom_ribbon(data = plot_frame, 
              aes(ymin = lower, ymax =upper), 
              alpha = 0.3) +
  ylab("Expected Count") +
  xlab("Week") +
  theme_light()

