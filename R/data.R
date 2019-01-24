#' Data frame of available eBird Status and Trends species
#'
#' A dataset containing the species for which eBird Status and Trends data are
#' avialable. In addition, the dates defining the boundaries of the seasons are
#' provided. These seasons are defined on a species-specific basis through
#' expert review. For information on the details of defining seasons, please see
#' the [seasons section of the
#' FAQ](https://ebird.org/science/status-and-trends/faq#seasons). Note that
#' missing dates imply that a season failed expert review for that species
#' within that season.
#'
#' @format A data frame with 107 rows and 14 variables:
#' \describe{
#'   \item{species_code}{Six letter eBird code in eBird Taxonomy v2016}
#'   \item{run_name}{Unique analysis identifier and the top level folder name
#'                   for all results}
#'   \item{scientific_name}{Scientific name from eBird Taxonomy v2016}
#'   \item{common_name}{English common name from eBird Taxonomy v2016}
#'   \item{breeding_start_dt}{Breeding season start date}
#'   \item{breeding_end_dt}{Breeding season start date}
#'   \item{nonbreeding_start_dt}{Non-breeding season start date}
#'   \item{nonbreeding_end_dt}{Non-breeding season start date}
#'   \item{postbreeding_migration_start_dt}{Post-breeding season start date}
#'   \item{postbreeding_migration_end_dt}{Post-breeding season start date}
#'   \item{prebreeding_migration_start_dt}{Pre-breeding season start date}
#'   \item{prebreeding_migration_end_dt}{Pre-breeding season start date}
#'   \item{year_round_start_dt}{For resident species, the year-round start date}
#'   \item{year_round_end_dt}{For resident species, the year-round end date}
#' }
"ebirdst_runs"

#' eBird Status and Trends predictors
#'
#' A data frame of the predictors used in the eBird Status and Trends models.
#' These include effort variables (e.g. distance travelled, number of observers,
#' etc.) in addition to land and water cover variables. These landcover
#' variables are derived from the MODIS MCD12Q1 500 m landcover product, and for
#' each land cover class four FRAGSTATS metrics are calculated within a 1.5 km
#' buffer around each checklist: % landcover (PLAND), edge density (ED), largest
#' patch index (LPI), and patch density (PD).
#'
#' @format A data frame with 87 rows and 5 columns:
#' \describe{
#'   \item{predictor}{Predictor variable name.}
#'   \item{predictor_tidy}{Predictor variable name, tidied to only contain
#'         lowercase letters and underscores.}
#'   \item{predictor_label}{Descriptive labels for predictors for plotting and
#'         translating the cryptic variables names (e.g. `umd_fs_c1` is
#'         Evergreen Needleleaf Forest.}
#'   \item{lc_class}{For the land and water cover FRAGSTATS variables, this
#'         gives the associated landcover class. It can be used for grouping
#'         and summarizing the four FRAGSTATS metrics to the level of the
#'         landcover class.}
#'   \item{lc_class_label}{Similar to `predictor_label`; however, this variable
#'         gives the four FRAGSTATS metrics a single name for the landcover
#'         class.}
#' }
"ebirdst_predictors"
