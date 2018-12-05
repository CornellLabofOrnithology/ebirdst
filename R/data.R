#' Data frame of species available
#'
#' A dataset containing the available species, with SPECIES_CODE, RUN_NAME,
#' SCI_NAME, and PRIMARY_COM_NAME.
#'
#' @format A data frame with 107 rows and 4 variables:
#' \describe{
#'   \item{SPECIES_CODE}{six letter eBird code in eBird Taxonomy v2016}
#'   \item{RUN_NAME}{unique analysis identifier and the top level folder name for all results}
#'   \item{SCI_NAME}{Scientific name from eBird Taxonomy v2016}
#'   \item{PRIMARY_COM_NAME}{English common name from eBird Taxonomy v2016}
#'   ...
#' }
"runs_w_names"
