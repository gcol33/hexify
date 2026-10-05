# Extracted from test-radius.R:193

# prequel ----------------------------------------------------------------------
MARS_RADIUS_KM <- 3389.50
MARS_SURFACE_KM2 <- 4 * pi * MARS_RADIUS_KM^2

# test -------------------------------------------------------------------------
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
hd <- hexify(df, lon = "lon", lat = "lat", resolution = 5, radius_km = "mars")
expect_equal(hexify:::grid_radius_km(hd@grid), MARS_RADIUS_KM)
expect_equal(as.data.frame(hd)$cell_area_km2,
               rep(hd@grid@area_km2, 3))
