# Clip hexagon grid to polygon boundary

Creates hexagon polygons clipped to a given polygon boundary. This is
useful for generating grids that conform to country borders, study
areas, or other irregular boundaries.

## Usage

``` r
grid_clip(boundary, grid, crop = TRUE)
```

## Arguments

- boundary:

  An sf/sfc polygon to clip to. Can be a country boundary, study area,
  or any polygon geometry.

- grid:

  A HexGridInfo object specifying the grid parameters

- crop:

  If TRUE (default), cells are cropped to the boundary. If FALSE, only
  cells whose centroids fall within the boundary are kept (no cropping).

## Value

sf object with hexagon polygons clipped to the boundary

## Details

The function first generates a rectangular grid covering the bounding
box of the input polygon, then clips or filters cells to the boundary.

When `crop = TRUE`, hexagons are geometrically intersected with the
boundary, which may produce partial hexagons at the edges. When
`crop = FALSE`, only complete hexagons whose centroids fall within the
boundary are returned.

## See also

[`grid_rect`](https://gcol33.github.io/hexify/reference/grid_rect.md)
for rectangular grids,
[`grid_global`](https://gcol33.github.io/hexify/reference/grid_global.md)
for global grids

## Examples

``` r
# \donttest{
# Get France boundary from built-in world map
france <- hexify_world[hexify_world$name == "France", ]

# Create grid clipped to France
grid <- hex_grid(area_km2 = 2000)
france_grid <- grid_clip(france, grid)

# Plot result
library(ggplot2)
ggplot() +
  geom_sf(data = france, fill = "gray95") +
  geom_sf(data = france_grid, fill = alpha("steelblue", 0.3),
          color = "steelblue") +
  theme_minimal()

# Keep only complete hexagons (no cropping)
france_grid_complete <- grid_clip(france, grid, crop = FALSE)
# }
```
