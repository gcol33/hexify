#' Simplified World Map
#'
#' A lightweight sf object containing simplified world country borders,
#' suitable for use as a basemap when visualizing hexagonal grids.
#'
#' @format An sf object with 177 features and 15 fields:
#' \describe{
#'   \item{name}{Country short name}
#'   \item{name_long}{Country full name}
#'   \item{admin}{Administrative name}
#'   \item{sovereignt}{Sovereignty}
#'   \item{iso_a2}{ISO 3166-1 alpha-2 country code}
#'   \item{iso_a3}{ISO 3166-1 alpha-3 country code}
#'   \item{iso_n3}{ISO 3166-1 numeric code}
#'   \item{continent}{Continent name}
#'   \item{region_un}{UN region}
#'   \item{subregion}{UN subregion}
#'   \item{region_wb}{World Bank region}
#'   \item{pop_est}{Population estimate}
#'   \item{gdp_md}{GDP in millions USD}
#'   \item{income_grp}{Income group classification}
#'   \item{economy}{Economy type}
#'   \item{geometry}{MULTIPOLYGON geometry in 'WGS84' (EPSG:4326)}
#' }
#'
#' @source Simplified from Natural Earth 1:110m Cultural Vectors
#'   (\url{https://www.naturalearthdata.com/})
#'
#' @examples
#' library(sf)
#'
#' # Plot the built-in world map
#' plot(st_geometry(hexify_world), col = "lightgray", border = "white")
#'
#' # Filter by continent
#' europe <- hexify_world[hexify_world$continent == "Europe", ]
#' plot(st_geometry(europe))
"hexify_world"
