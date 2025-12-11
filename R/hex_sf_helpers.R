# ---- Convert hex corners to sf polygon ----

#' Build an sf POLYGON from six (lon, lat) corner pairs
#'
#' @param lon numeric vector of length 6 (longitude)
#' @param lat numeric vector of length 6 (latitude)
#' @param crs integer CRS (default 4326)
#' @return sf object with one POLYGON geometry
#' @export
hex_corners_to_sf <- function(lon, lat, crs = 4326) {
  stopifnot(length(lon) == 6L, length(lat) == 6L)
  pts <- cbind(lon, lat)
  ring <- rbind(pts, pts[1, , drop = FALSE])  # close polygon
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = crs))
}
