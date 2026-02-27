# Run this script to regenerate test cache
# Usage: Rscript generate_cache.R

library(hexify)

cache <- list()

# Basic hexify results
df <- data.frame(lon = c(0, 10, 16.37, 2.35, -3.70), lat = c(45, 50, 48.21, 48.86, 40.42))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
cache$hexify_basic <- list(
  cell_id = result@cell_id,
  cell_cen_lon = result@cell_center[, "lon"],
  cell_cen_lat = result@cell_center[, "lat"],
  resolution = result@grid@resolution,
  aperture = result@grid@aperture
)

# Grid info
grid <- hex_grid(area_km2 = 1000)
cache$grid_1000 <- list(
  resolution = grid@resolution,
  aperture = grid@aperture,
  area_km2 = grid@area_km2,
  diagonal_km = grid@diagonal_km
)

# Cell conversions for aperture 3
cells3 <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 5, aperture = 3)
coords3 <- hexify_cell_to_lonlat(cells3, resolution = 5, aperture = 3)
cache$cells_ap3 <- list(cells = cells3, lon = coords3$lon_deg, lat = coords3$lat_deg)

# Cell conversions for aperture 4
cells4 <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 5, aperture = 4)
coords4 <- hexify_cell_to_lonlat(cells4, resolution = 5, aperture = 4)
cache$cells_ap4 <- list(cells = cells4, lon = coords4$lon_deg, lat = coords4$lat_deg)

# Cell conversions for aperture 7
cells7 <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 3, aperture = 7)
coords7 <- hexify_cell_to_lonlat(cells7, resolution = 3, aperture = 7)
cache$cells_ap7 <- list(cells = cells7, lon = coords7$lon_deg, lat = coords7$lat_deg)

# Projection forward
proj <- cpp_snyder_forward(10, 45)
cache$proj_forward <- as.list(proj)

# Face centers
cache$face_centers <- hexify_face_centers()

# Index conversions z3
cache$index_z3_sample <- hexify_cell_to_index(0, 5, 4, 3, 3, "z3")

# Quad IJ round-trip
qij <- hexify_cell_to_quad_ij(c(100, 1000, 5000), resolution = 10, aperture = 3)
cache$quad_ij <- qij

# Save
saveRDS(cache, "tests/testthat/fixtures/test_cache.rds")
cat("Cache saved with", length(cache), "entries\n")
