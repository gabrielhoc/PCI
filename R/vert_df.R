#' Threat data for 34124 global terrestrial vertebrates.
#'
#' A data set containing data on threat correlates for 34124 global terrestrial vertebrates.
#'
#' @format A data frame with 34124 rows and 18 variables:
#' \describe{
#'   \item{binomial}{character. species name}
#'   \item{class}{character. species class}
#'   \item{clim_2050_245}{numeric. Uninhabitable area due to climate change in 2050, under ssp 245}
#'   \item{clim_2100_245}{numeric. Uninhabitable area due to climate change in 2100, under ssp 245}
#'   \item{clim_2050_585}{numeric. Uninhabitable area due to climate change in 2050, under ssp 585}
#'   \item{clim_2100_585}{numeric. Uninhabitable area due to climate change in 2100, under ssp 585}
#'   \item{landuse_2050_245}{numeric. Uninhabitable area due to land use change in 2050, under ssp 245}
#'   \item{landuse_2100_245}{numeric. Uninhabitable area due to land use change in 2100, under ssp 245}
#'   \item{landuse_2050_585}{numeric. Uninhabitable area due to land use change in 2050, under ssp 585}
#'   \item{landuse_2100_585}{numeric. Uninhabitable area due to land use change in 2100, under ssp 585}
#'   \item{range_area}{numeric. Protected area, in km2}
#'   \item{popdens_2050_245}{numeric. Mean population density in 2050, under ssp 245}
#'   \item{popdens_2100_245}{numeric. Mean population density in 2100, under ssp 245}
#'   \item{popdens_2050_585}{numeric. Mean population density in 2050, under ssp 585}
#'   \item{popdens_2100_585}{numeric. Mean population density in 2100, under ssp 585}
#'   \item{range_area}{numeric. Area of distribution range, in km2}
#'   \item{body_mass}{numeric. Maximum body mass, in grams}
#'   \item{iucn_cat}{character. IUCN Red List threat category in July 2022}
#'   ...
#' }
"vert_df"