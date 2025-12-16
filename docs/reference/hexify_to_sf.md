# Convert hexify result to sf object

Converts the output of
[`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md) to an
sf spatial object. Can return either point geometries (cell centers) or
polygon geometries (cell boundaries).

## Usage

``` r
hexify_to_sf(data, geometry = c("point", "polygon"), aperture = 3L, crs = 4326)
```

## Arguments

- data:

  Data frame returned by
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)
  containing cell_id, cell_cen_lon, cell_cen_lat columns

- geometry:

  Type of geometry: "point" for cell centers (default), "polygon" for
  cell boundaries

- aperture:

  Grid aperture: 3, 4, or 7. Only needed for polygon geometry.

- crs:

  Coordinate reference system (default 4326 = WGS84)

## Value

sf object with:

- All original columns from data

- geometry column (sfc_POINT or sfc_POLYGON)

## Details

This is the recommended way to convert hexify output to sf for spatial
operations or plotting with ggplot2/tmap.

For `geometry = "point"`: Creates point geometries at cell centers using
cell_cen_lon and cell_cen_lat. Fast, suitable for large datasets.

For `geometry = "polygon"`: Creates polygon geometries for cell
boundaries. Slower but useful for choropleth maps. Duplicate cells are
automatically handled (each cell boundary appears once).

## See also

[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) for the
main function,
[`hexify_to_polygons`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md)
for polygon-only conversion

Other sf conversion:
[`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md),
[`hexify_cell_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md),
[`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md),
[`hexify_to_polygons()`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)
library(sf)

# Create sample data
df <- data.frame(
  site = c("Paris", "Vienna", "Madrid"),
  lon = c(2.35, 16.37, -3.70),
  lat = c(48.86, 48.21, 40.42),
  value = c(100, 200, 150)
)

# Hexify the data
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# Convert to sf points (cell centers)
sf_points <- hexify_to_sf(result)
plot(sf_points["value"])

# Convert to sf polygons (cell boundaries)
sf_polys <- hexify_to_sf(result, geometry = "polygon")
plot(sf_polys["value"])

# Use with ggplot2
library(ggplot2)
ggplot(sf_polys) +
  geom_sf(aes(fill = value)) +
  theme_minimal()
} # }
```
