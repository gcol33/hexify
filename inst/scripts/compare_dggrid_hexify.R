# inst/scripts/compare_dggrid_hexify.R
suppressPackageStartupMessages({
  library(hexify)
  library(dggridR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---- config ---------------------------------------------------------------

# Sample lon/lat points (well inside faces; EU cities)
pts <- tibble(
  id  = paste0("p", 1:10),
  lon = c(16.3738,  6.1296, 14.4378, 11.5820,  8.5417,
          2.3522, 12.4964, -3.7038, 12.5683, 21.0122),
  lat = c(48.2082, 49.6116, 50.0755, 48.1351, 47.3769,
          48.8566, 41.9028, 40.4168, 55.6761, 52.2297)
)

areas_km2 <- c(863, 2591, 7774, 23322)
WRITE_CSV  <- FALSE

# ---- helpers --------------------------------------------------------------

# Haversine distance in meters (inputs in degrees)
hav_m <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi/180
  φ1 <- lat1 * to_rad; φ2 <- lat2 * to_rad
  dφ <- (lat2 - lat1) * to_rad
  dλ <- (lon2 - lon1) * to_rad
  a <- sin(dφ/2)^2 + cos(φ1) * cos(φ2) * sin(dλ/2)^2
  2 * 6371008.8 * asin(pmin(1, sqrt(a)))  # Earth mean radius
}

# hexify helper: Snyder → quantize Z3 at eff_res → hex center → inverse
# IMPORTANT: flip_classes = TRUE to match ISEA3H (dggridR).
hexify_center_ll <- function(lon_deg, lat_deg, eff_res) {
  fwd  <- cpp_snyder_forward(lon_deg, lat_deg)
  face <- as.integer(fwd[["face"]])
  tx   <- as.numeric(fwd[["tx"]])
  ty   <- as.numeric(fwd[["ty"]])

  digs <- cpp_hex_index_z3_quantize_digits(tx, ty, eff_res,
                                           center_thr = 0.4,
                                           flip_classes = TRUE)$digits
  cen  <- cpp_hex_index_z3_center(digs, flip_classes = TRUE)
  ll   <- hex_face_xy_to_ll(cen[["cx"]], cen[["cy"]], face)

  c(lon = as.numeric(ll[["lon"]]), lat = as.numeric(ll[["lat"]]))
}

# optional: tighten inverse precision (if exposed)
if ("hex_snyder_set_precision" %in% getNamespaceExports("hexify")) {
  hex_snyder_set_precision("ultra")
}

# ---- compare across areas -------------------------------------------------

all_comparisons <- vector("list", length(areas_km2))
names(all_comparisons) <- as.character(areas_km2)

for (area in areas_km2) {
  grid <- dgconstruct(area = area, topology = "HEXAGON", metric = TRUE, resround = "nearest")
  eff_res <- grid$res

  cat(sprintf(
    "Resolution: %d, Area (km^2): %.12f, Spacing (km): %.10f, CLS (km): %.10f\n",
    eff_res, grid$Area, grid$Spacing, grid$CLS
  ))

  # dggridR centers (coerce list -> tibble)
  seqnum <- dgGEO_to_SEQNUM(grid, pts$lon, pts$lat)$seqnum
  centers_list <- dgSEQNUM_to_GEO(grid, seqnum)  # list with lon_deg, lat_deg
  dgg_df <- tibble(
    id      = pts$id,
    dgg_lon = as.numeric(centers_list$lon_deg),
    dgg_lat = as.numeric(centers_list$lat_deg)
  )

  # hexify centers at same effective resolution
  hex_mat <- t(mapply(function(lo, la) hexify_center_ll(lo, la, eff_res), pts$lon, pts$lat))
  hex_df <- as_tibble(hex_mat) |>
    mutate(id = pts$id, .before = 1L) |>
    rename(hex_lon = lon, hex_lat = lat)

  # join and compute spherical distance error (meters)
  cmp <- dgg_df |>
    left_join(hex_df, by = "id") |>
    mutate(
      area_km2 = area,
      gc_err_m = hav_m(dgg_lon, dgg_lat, hex_lon, hex_lat)
    )

  all_comparisons[[as.character(area)]] <- cmp
}

cmp_all <- bind_rows(all_comparisons)

summary_tbl <- cmp_all |>
  group_by(area_km2) |>
  summarize(
    n           = n(),
    max_err_m   = max(gc_err_m, na.rm = TRUE),
    mean_err_m  = mean(gc_err_m, na.rm = TRUE),
    median_err_m= median(gc_err_m, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nGreat-circle error summary (hexify vs dggridR centers):\n")
print(summary_tbl, n = Inf)

if (isTRUE(WRITE_CSV)) {
  write_csv(cmp_all,    "compare_hexify_dggrid_points.csv")
  write_csv(summary_tbl,"compare_hexify_dggrid_summary.csv")
  cat("\nWrote: compare_hexify_dggrid_points.csv and compare_hexify_dggrid_summary.csv\n")
}
