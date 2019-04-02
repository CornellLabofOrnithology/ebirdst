#' Spatial grid sampling methods
#'
#' Methods for subsampling data to deal with spatiotemporal bias in
#' observations. These functions define a grid in space and time, then sample
#' the given number of points from each cell. [sample_case_control()]
#' additionally samples presence and absence independently.
#'
#' @param x data frame or [sf] point object; the points to subsample
#' @param res numeric; the size in meters of the grid to sample from. This can
#'   be a 2 element vector indicating the x and y dimensions of the cells.
#' @param t_res numeric; the temporal resolution for sampling expressed as a
#'   proportion of the year. For example, `7 / 375` would result in sampling
#'   from each week.
#' @param n integer; the number of points to sample from each grid cell.
#' @param replace logical; whether to sample with replacement.
#' @param jitter logical; to avoid always using the same grid for sampling,
#'   the grid can be jittered so that the origin is different each time this
#'   function is called.
#'
#' @return Logical vector indicating which rows are selected.
#'
#' @export
#' @rdname ebirdst_sample
#'
#' @examples
#' # download example data
#' sp_path <- ebirdst_download("example_data", tifs_only = FALSE)
#'
#' # test data to sample
#' test_data <- load_test_data(sp_path, return_sf = TRUE)
#'
#' # sample on a 100km, 1 month grid
#' s <- sample_grid(test_data, res = 100000, t_res = 1 / 12)
#' td_grid <- test_data[s, ]
#'
#' # case control sampling independently samples presence and absence
#' s <- sample_case_control(test_data, res = 100000, t_res = 1 / 12)
#' td_cc <- test_data[s, ]
#'
#' # grid sampling preserves the presence/absence ratio
#' table(test_data$obs > 0) / nrow(test_data)
#' table(td_grid$obs > 0) / nrow(td_grid)
#' # while case control sampling increases the prevelance of presences
#' table(td_cc$obs > 0) / nrow(td_cc)
#'
#' \dontrun{
#'
#' # plot
#' library(sf)
#' p <- par(mar = c(0, 0, 0, 0))
#' plot(st_geometry(test_data), col = "black", pch = 19, cex = 0.2)
#' plot(st_geometry(td_cc), col = "red", pch = 19, cex = 0.5, add = TRUE)
#' par(p)
#'
#' }
sample_grid <- function(x, res, t_res, n = 1, replace = FALSE,
                        jitter = TRUE) {
  UseMethod("sample_grid")
}

#' @export
sample_grid.sf <- function(x, res, t_res, n = 1, replace = FALSE,
                           jitter = TRUE) {
  stopifnot(nrow(x) > 0, "date" %in% names(x),
            sf::st_geometry_type(x) == "POINT")
  stopifnot(is.numeric(res), length(res) %in% 1:2,
            all(!is.na(res)), all(res > 1))
  stopifnot(is.numeric(t_res), length(t_res) == 1, !is.na(res),
            t_res >= 0, t_res <= 1)
  stopifnot(is_integer(n), length(n) == 1, !is.na(n), n > 0)
  stopifnot(is.logical(replace), length(replace) == 1, !is.na(replace))
  stopifnot(is.logical(jitter), length(jitter) == 1, !is.na(jitter))

  if (length(res) == 1) {
    res <- rep(res, 2)
  }

  # use sinusoidal equal area projection
  sinu <- paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                "+a=6371007.181 +b=6371007.181 +units=m +no_defs")
  x <- sf::st_transform(x, crs = sinu)
  x <- cbind(sf::st_coordinates(x), x$date)
  x <- stats::setNames(as.data.frame(x), c("x", "y", "date"))

  # define grid
  # lower left corner
  ll <- apply(x, 2, min, na.rm = TRUE)
  # jitter
  if (jitter) {
    ll <- ll - stats::runif(3) * c(res, t_res)
  }

  # assign to grid cells
  x_cell <- 1 + (x$x - ll[1]) %/% res[1]
  y_cell <- 1 + (x$y - ll[2]) %/% res[2]
  t_cell <- 1 + (x$date - ll[3]) %/% t_res
  cell <- x_cell +
    (y_cell - 1) * max(x_cell, na.rm = TRUE) +
    (t_cell - 1) * max(x_cell, na.rm = TRUE) * max(y_cell, na.rm = TRUE)
  cell <- as.factor(as.integer(cell))

  # sample from grid cells
  sampled <- tapply(seq_along(cell)[!is.na(cell)],
                    cell[!is.na(cell)],
                    safe_sample, size = n, replace = replace,
                    simplify = FALSE)
  sampled <- do.call(c, sampled)
  is_sampled <- rep(FALSE, nrow(x))
  is_sampled[sampled] <- TRUE

  return(is_sampled)
}

#' @export
sample_grid.data.frame <- function(x, res, t_res, n = 1, replace = FALSE,
                                   jitter = TRUE) {
  stopifnot(nrow(x) > 0, all(c("lon", "lat", "date") %in% names(x)))
  x <- sf::st_as_sf(x, coords = c("lon", "lat"), crs = 4326)
  sample_grid.sf(x, res = res, t_res = t_res, n = n,
                 replace = replace, jitter = jitter)
}


#' @export
#' @rdname ebirdst_sample
sample_case_control <- function(x, res, t_res, n = 1, replace = FALSE,
                                jitter = TRUE) {
  UseMethod("sample_case_control")
}


#' @export
sample_case_control.sf <- function(x, res, t_res, n = 1, replace = FALSE,
                                   jitter = TRUE)  {
  stopifnot(nrow(x) > 1, all(c("obs", "date") %in% names(x)),
            sf::st_geometry_type(x) == "POINT")
  stopifnot(is.numeric(res), length(res) %in% 1:2,
            all(!is.na(res)), all(res > 1))
  stopifnot(is.numeric(t_res), length(t_res) == 1, !is.na(res),
            t_res >= 0, t_res <= 1)
  stopifnot(is_integer(n), length(n) == 1, !is.na(n), n > 0)
  stopifnot(is.logical(replace), length(replace) == 1, !is.na(replace))
  stopifnot(is.logical(jitter), length(jitter) == 1, !is.na(jitter))

  # spit into presence and absence
  x$row_id <- seq_len(nrow(x))
  x_split <- split(x, ifelse(x$obs > 0, "pres", "abs"))

  # sample indedendently
  for (i_pa in names(x_split)) {
    to_keep <- sample_grid(x_split[[i_pa]], res = res, t_res = t_res, n = n,
                           replace = replace, jitter = jitter)
    x_split[[i_pa]] <- x_split[[i_pa]]$row_id[to_keep]
  }

  # combine
  sampled <- do.call(c, x_split)
  is_sampled <- rep(FALSE, nrow(x))
  is_sampled[sampled] <- TRUE

  return(is_sampled)
}


#' @export
sample_case_control.data.frame <- function(x, res, t_res, n = 1,
                                           replace = FALSE, jitter = TRUE) {
  stopifnot(nrow(x) > 0, all(c("obs", "lon", "lat", "date") %in% names(x)))
  x <- sf::st_as_sf(x, coords = c("lon", "lat"), crs = 4326)
  sample_case_control.sf(x, res = res, t_res = t_res, n = n,
                         replace = replace, jitter = jitter)
}

safe_sample <- function(x, size, ...) {
  stopifnot(is.numeric(size), length(size) == 1, size >= 1)
  if (length(x) <= size || length(x) == 1) {
    return(x)
  }
  sample(x, size = size, ...)
}
