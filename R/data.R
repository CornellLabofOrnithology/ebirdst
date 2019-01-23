#' Data frame of species available
#'
#' A dataset containing the available species, with SPECIES_CODE, RUN_NAME,
#' SCI_NAME, and PRIMARY_COM_NAME. In addition, the dates defining the
#' boundaries of the seasons are provided. These seasons are defined on a
#' species-specific basis through expert review. For information on the details
#' of defining seasons, please see the [seasons section of the
#' FAQ](https://ebird.org/science/status-and-trends/faq#seasons). Note that
#' missing dates imply that a season failed expert review for that species
#' within that season.
#'
#' @format A data frame with 107 rows and 14 variables:
#' \describe{
#'   \item{SPECIES_CODE}{Six letter eBird code in eBird Taxonomy v2016}
#'   \item{RUN_NAME}{Unique analysis identifier and the top level folder name
#'                   for all results}
#'   \item{SCI_NAME}{Scientific name from eBird Taxonomy v2016}
#'   \item{PRIMARY_COM_NAME}{English common name from eBird Taxonomy v2016}
#'   \item{BREEDING_START_DT}{Breeding season start date}
#'   \item{BREEDING_END_DT}{Breeding season start date}
#'   \item{NONBREEDING_START_DT}{Non-breeding season start date}
#'   \item{NONBREEDING_END_DT}{Non-breeding season start date}
#'   \item{POSTBREEDING_MIGRATION_START_DT}{Post-breeding season start date}
#'   \item{POSTBREEDING_MIGRATION_END_DT}{Post-breeding season start date}
#'   \item{PREBREEDING_MIGRATION_START_DT}{Pre-breeding season start date}
#'   \item{PREBREEDING_MIGRATION_END_DT}{Pre-breeding season start date}
#'   \item{YEAR_ROUND_START_DT}{For resident species, the year-round start date}
#'   \item{YEAR_ROUND_END_DT}{For resident species, the year-round end date}
#' }
"runs_w_names"

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
