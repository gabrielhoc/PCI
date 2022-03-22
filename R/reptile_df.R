#' Threat data for 10851 global reptiles.
#'
#' A dataset containing data on threat correlates for 10851 global reptiles
#'
#' @format A data frame with 359 rows and 3 variables:
#' \describe{
#'   \item{sp}{character. species name}
#'   \item{iucn_cat}{character. IUCN Red List threat category}
#'   \item{insular}{binary. If insular endemic, value 1, if not, value 0}
#'   \item{range_area}{numeric. Area of distribution range}
#'   \item{human_footprint}{numeric. Average human footprint (ref) density across species distribution}
#'   \item{pop_dens}{numeric. Average human population density across species distribution}
#'   ...
#' }
"reptile_df"
