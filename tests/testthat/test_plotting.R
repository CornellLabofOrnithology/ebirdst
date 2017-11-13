context("Plotting functions")

# plot_pis
test_that("stemhelper plot_pis", {
  root_path <- "~/Box Sync/Projects/2015_stem_hwf/documentation/data-raw/"
  species <- "woothr-ERD2016-PROD-20170505-3f880822"
  sp_path <- paste(root_path, species, sep = "")

  pis <- load_pis(sp_path)

  # expected with st_extent
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        by_cover_class = TRUE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 10),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = TRUE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        pretty_names = FALSE),
               NA)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = TRUE,
                        print_plot = FALSE),
               NA)

  # checking return quality
  expect_equal(length(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 25,
                        return_top = TRUE)),
               25)
  expect_is(plot_pis(path = sp_path,
                     pis = pis,
                     st_extent = ne_extent,
                     num_top_preds = 25,
                     return_top = TRUE),
               "character")

  # nothing to do
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        return_top = FALSE,
                        print_plot = FALSE),
               "Both print and return params are FALSE")

  # not enough num_top_preds
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent,
                        num_top_preds = 1),
               "num_top_preds must be greater than 1")

  # broken path
  sp_path <- '~/some/messed/up/path/that/does/not/exist'
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "RData file does not exist")

  sp_path <- paste(root_path, species, sep = "")

  # missing temporal info
  ne_extent <- list(type = "rectangle",
                    lat.min = 40,
                    lat.max = 47,
                    lon.min = -80,
                    lon.max = -70)
  expect_error(plot_pis(path = sp_path,
                     pis = pis,
                     st_extent = ne_extent),
               "Must provide t.min and t.max")

  # reversed min max
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    lon.max = -70,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "Latitude maximum is less than latitude minimum")

  # missing a corner
  ne_extent <- list(type = "rectangle",
                    lat.min = 47,
                    lat.max = 40,
                    lon.min = -80,
                    t.min = 0.425,
                    t.max = 0.475)
  expect_error(plot_pis(path = sp_path,
                        pis = pis,
                        st_extent = ne_extent),
               "Missing max longitude")

  # st_extent is not list
  ne_extent <- c(type = "rectangle",
                 lat.min = 40,
                 lat.max = 47,
                 lon.min = -80,
                 lon.max = -70,
                 t.min = 0.425,
                 t.max = 0.475)
  expect_error(stack_stem(sp_path,
                          variable = "abundance_umean",
                          st_extent = ne_extent),
               "st_extent argument must be a list")
})

# plot_pds
# cake_plots
