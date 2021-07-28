#' #' Calculate summary statistics over a set of regional polygons
#' #'
#' #' Given a set of polygons defining regions, calculate a suite of statistics
#' #' for the range and abundance within those regions.
#' #' @inheritParams load_raster
#' #' @param regions [sf] object;
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' \donttest{
#' #' library(sf)
#' #'
#' #' # download example data
#' #' path <- ebirdst_download("example_data", tifs_only = FALSE)
#' #' # or get the path if you already have the data downloaded
#' #' path <- get_species_path("example_data")
#' #'
#' #' # generate an example "region"
#' #' ex_region <- st_point(c(-84.6, 44.3)) %>%
#' #'   st_sfc(crs = 4326) %>%
#' #'   st_transform(crs = "+proj=sinu") %>%
#' #'   st_buffer(dist = 10000) %>%
#' #'   st_sf(region_name = "Buffered Point", geometry = .)
#' #'
#' #'
#' #' }
#' regional_stats <- function(path, regions) {
#'   stopifnot(is.character(path), length(path) == 1, dir.exists(path))
#'   stopifnot(inherits(regions, "sf"))
#'
#'   # load data
#'   p <- load_config(path)
#'   is_resident <- p$IS_RESIDENT
#'   abd_w <- load_raster(path, product = "abundance", res = "lr")
#'   abd_s <- load_raster(path, product = "abundance_seasonal", res = "lr")
#'   if (is_resident) {
#'     complete_seasons <- "resident"
#'   } else {
#'     complete_seasons <- c("nonbreeding", "prebreeding_migration",
#'                           "breeding", "postbreeding_migration")
#'   }
#'   complete_seasons <- intersect(names(abd_s), complete_seasons)
#'   if (length(complete_seasons) == 0) {
#'     stop("No valid seasonal abundance layers found")
#'   }
#'   abd_s <- abd_s[[complete_seasons]]
#'
#'   # assign unique id to regions
#'   regions <- st_transform(regions, crs = raster::projection(abd_s))
#'   regions_df <- sf::st_drop_geometry(regions)
#'   regions[[".region_id"]] <- seq_len(nrow(regions))
#'   regions <- regions[, ".region_id"]
#'
#'   # find regions that have predictions
#'   has_preds <- raster_summarize(abd_s, regions, fun = "count")
#'   has_preds$prediction <- has_preds$value > 0
#'   has_preds$value <- NULL
#'
#'   # use days
#'   # determine weekly pct region occupied
#'   region_occupied <- function(wk) {
#'     pa_wk <- r_week[[wk]] > 0
#'     pa_wk[is.na(pa_wk[])] <- 0
#'     exactextractr::exact_extract(abd_w, regions, fun = "mean") %>%
#'       tibble(v = .) %>%
#'       transmute(season = week_seasons[wk],
#'                 # is the region occupied for that week
#'                 occupied = coalesce(v, 0) >= threshold_occupied) %>%
#'       bind_cols(regions_df, .)
#'   }
#'   valid_weeks <- which(!is.na(week_seasons))
#'   region_occ <- mclapply_check(X = valid_weeks,
#'                                FUN = region_occupied,
#'                                mc.preschedule = TRUE, mc.set.seed = TRUE,
#'                                mc.cores = NUM_CORES) %>%
#'     bind_rows()
#'   collect_garbage()
#'   # use days
#'   ud_regions <- region_occ %>%
#'     filter(!is.na(season)) %>%
#'     group_by(region_type, region_code, season) %>%
#'     summarise(days_occ = sum(occupied), .groups = "drop") %>%
#'     mutate(variable = "range_days_occupation",
#'            value = round(7 * days_occ)) %>%
#'     select(region_type, region_code, season, variable, value)
#'
#' }
#'
#' raster_summarize <- function(x, y, fun) {
#'   e <- exactextractr::exact_extract(x = x, y = y, fun = fun,
#'                                     force_df = TRUE,
#'                                     full_colnames =TRUE,
#'                                     append_cols = TRUE,
#'                                     progress = FALSE)
#'   e <- tidyr::pivot_longer(data = e, cols = dplyr::starts_with(fun),
#'                            names_to = "season",
#'                            names_prefix = paste0(fun, "."),
#'                            values_to = "value")
#'   return(e)
#' }
#'
#' raster_extract <- function(x, y) {
#'   e <- exactextractr::exact_extract(x = x, y = y, fun = fun,
#'                                     force_df = TRUE,
#'                                     full_colnames =TRUE,
#'                                     append_cols = TRUE,
#'                                     progress = FALSE)
#'   e <- tidyr::pivot_longer(data = e, cols = dplyr::starts_with(fun),
#'                            names_to = "season",
#'                            names_prefix = paste0(fun, "."),
#'                            values_to = "value")
#'   return(e)
#' }
