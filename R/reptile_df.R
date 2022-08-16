#' Threat data for 10851 global reptiles.
#'
#' A dataset containing data on threat correlates for 10851 global reptiles
#'
#' @format A data frame with 359 rows and 3 variables:
#' \describe{
#'   \item{sp}{character. species name}
#'   \item{iucn_cat}{character. IUCN Red List threat category}
#'   \item{mass}{numeric. Maximum body mass, in grams}
#'   \item{brood_size}{numeric. Average brood size}
#'   \item{insular}{binary. If insular endemic, value 1, if not, value 0}
#'   \item{range_area}{numeric. Area of distribution range, in km2}
#'   \item{human_footprint}{numeric. Average human footprint (ref) density across species distribution}
#'   \item{clim_2050_245}{numeric. Uninhabitable area due to climate change in 2050, under ssp 245}
#'   \item{clim_2100_245}{numeric. Uninhabitable area due to climate change in 2100, under ssp 245}
#'   \item{clim_2050_585}{numeric. Uninhabitable area due to climate change in 2050, under ssp 585}
#'   \item{clim_2100_585}{numeric. Uninhabitable area due to climate change in 2100, under ssp 585}
#'   \item{luh_2050_245}{numeric. Uninhabitable area due to land use change in 2050, under ssp 245}
#'   \item{luh_2100_245}{numeric. Uninhabitable area due to land use change in 2100, under ssp 245}
#'   \item{luh_2050_585}{numeric. Uninhabitable area due to land use change in 2050, under ssp 585}
#'   \item{luh_2100_585}{numeric. Uninhabitable area due to land use change in 2100, under ssp 585}
#'   ...
#' }
"reptile_df"
