#' Assign hex cells (ISEA3H, aperture 3) for lon/lat
#'
#' Forward -> quantize (Z3) -> center (face) -> inverse to lon/lat.
#' Optionally return polygons (sf), in which case sf must be installed.
#'
#' @param lon,lat numeric vectors (same length), degrees.
#' @param eff_res integer effective resolution (>= 1).
#' @param flip_classes logical; TRUE matches ISEA3H parity used by dggridR.
#' @param make_polygons logical; if TRUE, return an sf with hex polygons.
#' @return data.frame with id, face, eff_res, center_lon, center_lat;
#'         if make_polygons=TRUE, an sf with geometry column.
#' @export
hex_assign <- function(lon, lat, eff_res,
                       flip_classes = TRUE,
                       make_polygons = FALSE) {
  stopifnot(length(lon) == length(lat), length(eff_res) == 1L, eff_res >= 1L)

  n <- length(lon)

  # Ensure ISEA3H default orientation (safe no-op if already set)
  if ("cpp_build_icosa" %in% getNamespaceExports("hexify")) {
    cpp_build_icosa(11.25, 58.28252559, 0)
  }
  hexify_set_precision("ultra")

  # 1) Forward
  fwd <- lapply(seq_len(n), function(k) cpp_snyder_forward(lon[k], lat[k]))
  face <- vapply(fwd, function(v) as.integer(v[["face"]]), integer(1))
  tx   <- vapply(fwd, function(v) as.numeric(v[["tx"]]),   numeric(1))
  ty   <- vapply(fwd, function(v) as.numeric(v[["ty"]]),   numeric(1))

  # 2) Quantize digits
  digs <- lapply(seq_len(n), function(k)
    cpp_hex_index_z3_quantize_digits(tx[k], ty[k], eff_res,
                                     center_thr = 0.4,
                                     flip_classes = flip_classes)$digits)

  # Canonical string id: Face + digits
  id <- vapply(seq_len(n), function(k) {
    paste0("F", face[k], ":Z3:", paste(digs[[k]], collapse = "-"))
  }, character(1))

  # 3) Center (face) -> lon/lat
  cen <- lapply(digs, function(d) cpp_hex_index_z3_center(d, flip_classes = flip_classes))
  cx  <- vapply(cen, function(c) as.numeric(c[["cx"]]), numeric(1))
  cy  <- vapply(cen, function(c) as.numeric(c[["cy"]]), numeric(1))
  llc <- lapply(seq_len(n), function(k) hexify_inverse(cx[k], cy[k], face[k]))
  clon <- vapply(llc, function(p) as.numeric(p[["lon"]]), numeric(1))
  clat <- vapply(llc, function(p) as.numeric(p[["lat"]]), numeric(1))

  df <- data.frame(
    id, face, eff_res,
    center_lon = clon,
    center_lat = clat,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(make_polygons)) return(df)

  # 4) Polygons (sf)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("make_polygons=TRUE requires the 'sf' package.")
  }
  polys <- vector("list", n)
  for (k in seq_len(n)) {
    cor <- cpp_hex_index_z3_corners(digs[[k]],
                                    flip_classes = flip_classes,
                                    hex_radius = 1.0)
    lls <- vapply(1:6, function(i) hexify_inverse(cor$x[i], cor$y[i], face[k]),
                  c(lon = 0.0, lat = 0.0))
    # build a small sf for each hex
    polys[[k]] <- hex_corners_to_sf(lls["lon", ], lls["lat", ])
  }

  # extract sfgs and build a single sfc with CRS
  sfg_list <- lapply(polys, function(p) sf::st_geometry(p)[[1]])
  sfc      <- sf::st_sfc(sfg_list, crs = 4326)

  sf::st_sf(df, geometry = sfc)

}
