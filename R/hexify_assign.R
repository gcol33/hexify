# =============================================================================
# POINT ASSIGNMENT (ISEA3H, APERTURE 3)
# =============================================================================

#' Assign hex cells ('ISEA3H', aperture 3) for lon/lat
#'
#' Forward projection -> quantize -> cell centre -> inverse to lon/lat, on the
#' aperture-3 pipeline used by [hexify()]: points fold into the non-negative
#' quad frame ([lonlat_to_cell()]), and the centre comes back from the cell
#' itself ([cell_to_lonlat()]), so the returned centre is the centre of the
#' assigned cell.
#'
#' @param lon,lat numeric vectors (same length), degrees.
#' @param effective_res integer effective resolution (>= 1).
#' @param make_polygons logical; if TRUE, return an sf with hex polygons.
#' @return data.frame with id (Z3 index string), face (quad, 0-11),
#'         effective_res, center_lon, center_lat; if make_polygons=TRUE,
#'         an sf with geometry column.
#'
#' @keywords internal
#' @export
hexify_assign <- function(lon, lat, effective_res, make_polygons = FALSE) {

  stopifnot(
    length(lon) == length(lat),
    length(effective_res) == 1L,
    effective_res >= 1L
  )

  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  res <- as.integer(effective_res)

  cell_id <- cpp_lonlat_to_cell(lon, lat, res, 3L)
  center <- cpp_cell_to_lonlat(cell_id, res, 3L)
  quad_ij <- cpp_cell_to_quad_ij(cell_id, res, 3L)

  df <- data.frame(
    id = vapply(seq_along(cell_id), function(k) {
      cpp_cell_to_index(quad_ij$quad[k], quad_ij$i[k], quad_ij$j[k],
                        res, 3L, "z3")
    }, character(1)),
    face = quad_ij$quad,
    effective_res = res,
    center_lon = center$lon_deg,
    center_lat = center$lat_deg,
    stringsAsFactors = FALSE
  )

  if (!isTRUE(make_polygons)) return(df)

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("make_polygons=TRUE requires the 'sf' package.")
  }

  sf::st_sf(df, geometry = isea_cells_to_sfc(cell_id, res, 3L))
}
