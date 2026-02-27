# =============================================================================
# INTERNAL HELPERS FOR HEX ASSIGNMENT
# =============================================================================

#' Project points and quantize to cell structure
#' @noRd
project_and_quantize <- function(lon, lat, effective_res, match_dggrid_parity) {
  n <- length(lon)

  # Forward projection: lon/lat -> (face, icosa_triangle_x, icosa_triangle_y)
  proj <- cpp_snyder_forward(lon, lat)
  face <- as.integer(proj[["face"]])
  icosa_triangle_x <- as.numeric(proj[["icosa_triangle_x"]])
  icosa_triangle_y <- as.numeric(proj[["icosa_triangle_y"]])

  # Quantize to cell digits
  quant <- cpp_hex_index_z3_quantize_digits(
    icosa_triangle_x, icosa_triangle_y, effective_res,
    center_thr = 0.4,
    flip_classes = match_dggrid_parity
  )

  list(
    face = face,
    icosa_triangle_x = icosa_triangle_x,
    icosa_triangle_y = icosa_triangle_y,
    digits = quant$digits
  )
}

#' Compute cell center coordinates from digits
#' @noRd
compute_cell_center <- function(digits, face, match_dggrid_parity) {
  center <- cpp_hex_index_z3_center(digits, flip_classes = match_dggrid_parity)
  coords <- hexify_inverse(center[["cx"]], center[["cy"]], face)
  c(lon = coords[["lon"]], lat = coords[["lat"]])
}

#' Generate canonical cell ID string
#' @noRd
make_cell_id <- function(face, digits) {
  paste0("F", face, ":Z3:", paste(digits, collapse = "-"))
}

#' Build hex polygon from cell digits
#' @noRd
build_hex_polygon <- function(digits, face, match_dggrid_parity) {
  hex_corners <- cpp_hex_index_z3_corners(
    digits,
    flip_classes = match_dggrid_parity,
    hex_radius = 1.0
  )

  corner_lonlat <- vapply(1:6, function(i) {
    hexify_inverse(hex_corners$x[i], hex_corners$y[i], face)
  }, c(lon = 0.0, lat = 0.0))

  hex_corners_to_sf(corner_lonlat["lon", ], corner_lonlat["lat", ])
}

# =============================================================================
# PUBLIC API
# =============================================================================

#' Assign hex cells ('ISEA3H', aperture 3) for lon/lat
#'
#' Forward -> quantize (Z3) -> center (face) -> inverse to lon/lat.
#' Optionally return polygons (sf), in which case sf must be installed.
#'
#' @param lon,lat numeric vectors (same length), degrees.
#' @param effective_res integer effective resolution (>= 1).
#' @param match_dggrid_parity logical; TRUE matches 'ISEA3H' parity used by 'dggridR'.
#' @param make_polygons logical; if TRUE, return an sf with hex polygons.
#' @return data.frame with id, face, effective_res, center_lon, center_lat;
#'         if make_polygons=TRUE, an sf with geometry column.
#'
#' @keywords internal
#' @export
hexify_assign <- function(lon, lat, effective_res,
                          match_dggrid_parity = TRUE,
                          make_polygons = FALSE) {

  stopifnot(
    length(lon) == length(lat),
    length(effective_res) == 1L,
    effective_res >= 1L
  )

  n <- length(lon)

  # Initialize icosahedron with ISEA3H default orientation
  if ("cpp_build_icosa" %in% getNamespaceExports("hexify")) {
    cpp_build_icosa(ISEA_VERT0_LON_DEG, ISEA_VERT0_LAT_DEG, ISEA_AZIMUTH_DEG)
  }
  hexify_set_precision("ultra")

  # Process each point: project, quantize, compute centers
  results <- lapply(seq_len(n), function(k) {
    proj <- project_and_quantize(
      lon[k], lat[k], effective_res, match_dggrid_parity
    )
    center <- compute_cell_center(
      proj$digits, proj$face, match_dggrid_parity
    )

    list(
      id = make_cell_id(proj$face, proj$digits),
      face = proj$face,
      digits = proj$digits,
      center_lon = center["lon"],
      center_lat = center["lat"]
    )
  })

  # Extract vectors from results
  df <- data.frame(
    id = vapply(results, `[[`, character(1), "id"),
    face = vapply(results, `[[`, integer(1), "face"),
    effective_res = effective_res,
    center_lon = vapply(results, `[[`, numeric(1), "center_lon"),
    center_lat = vapply(results, `[[`, numeric(1), "center_lat")
  )

  if (!isTRUE(make_polygons)) return(df)

  # Build polygons (requires sf)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("make_polygons=TRUE requires the 'sf' package.")
  }

  sfg_list <- lapply(seq_len(n), function(k) {
    poly <- build_hex_polygon(
      results[[k]]$digits, results[[k]]$face, match_dggrid_parity
    )
    sf::st_geometry(poly)[[1]]
  })

  sfc <- sf::st_sfc(sfg_list, crs = 4326)
  sf::st_sf(df, geometry = sfc)
}
