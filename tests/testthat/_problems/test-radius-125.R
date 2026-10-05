# Extracted from test-radius.R:125

# prequel ----------------------------------------------------------------------
MARS_RADIUS_KM <- 3389.50
MARS_SURFACE_KM2 <- 4 * pi * MARS_RADIUS_KM^2

# test -------------------------------------------------------------------------
mars <- hex_grid(resolution = 5, radius_km = "mars")
cells <- lonlat_to_cell(c(0, 30), c(0, 40), mars)
earth <- hex_grid(resolution = 5)
expect_equal(unname(cell_area(cells, mars)),
               unname(cell_area(cells, earth)) * (MARS_RADIUS_KM / EARTH_RADIUS_KM)^2)
