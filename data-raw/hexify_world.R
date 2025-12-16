# data-raw/hexify_world.R
# Script to generate the built-in simplified world map
#
# This creates a lightweight world borders sf object from Natural Earth data.
# The map is simplified to keep the package small.
#
# To regenerate: source("data-raw/hexify_world.R")

library(sf)
library(rnaturalearth)

# Disable S2 for simplification (avoids edge crossing errors)
sf_use_s2(FALSE)

# Download Natural Earth 50m scale (medium resolution - good balance)
world <- ne_countries(scale = 50, returnclass = "sf")

# Keep useful columns (identification, geography, economy)
keep_cols <- c(
  # Names
  "name",           # Short name
  "name_long",      # Full name
  "admin",          # Administrative name
  "sovereignt",     # Sovereignty
  # ISO codes
  "iso_a2",         # ISO 3166-1 alpha-2
  "iso_a3",         # ISO 3166-1 alpha-3
  "iso_n3",         # ISO 3166-1 numeric
  # Geography
  "continent",
  "region_un",      # UN region
  "subregion",      # UN subregion
  "region_wb",      # World Bank region
  # Economy/population
  "pop_est",        # Population estimate
  "gdp_md",         # GDP in millions USD
  "income_grp",     # Income group
  "economy"         # Economy type
)
world <- world[, keep_cols]

# Simplify geometry (tolerance in degrees, ~10km at equator)
world <- st_simplify(world, dTolerance = 0.1, preserveTopology = TRUE)

# Ensure valid geometries
world <- st_make_valid(world)

# Keep only land polygons (remove any empty geometries)
world <- world[!st_is_empty(world), ]

# Set standard CRS
st_crs(world) <- 4326

# Rename to hexify_world
hexify_world <- world

# Save as internal data
usethis::use_data(hexify_world, overwrite = TRUE)

# Report size
cat("hexify_world size:", format(object.size(hexify_world), units = "KB"), "\n")
cat("Number of features:", nrow(hexify_world), "\n")
cat("Columns:", paste(names(hexify_world), collapse = ", "), "\n")
