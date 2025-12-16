# Generate a rectangular grid of hexagon polygons

Creates hexagon polygons covering a rectangular geographic region.

## Usage

``` r
hexify_grid_rect(
  minlon,
  maxlon,
  minlat,
  maxlat,
  area,
  aperture = 3L,
  resround = "nearest"
)
```

## Arguments

- minlon, maxlon:

  Longitude bounds

- minlat, maxlat:

  Latitude bounds

- area:

  Target cell area in km²

- aperture:

  Grid aperture: 3, 4, or 7

- resround:

  Resolution rounding: "nearest", "up", or "down"

## Value

sf object with hexagon polygons covering the specified region

## See also

[`hexify_grid_global`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md)
for global grids,
[`hexify_to_polygons`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md)
for data-driven polygon generation

Other sf conversion:
[`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md),
[`hexify_cell_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md),
[`hexify_to_polygons()`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md),
[`hexify_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)
library(sf)

grid <- hexify_grid_rect(
  minlon = -10, maxlon = 20,
  minlat = 35, maxlat = 60,
  area = 5000
)
plot(st_geometry(grid), border = "gray")
} # }
```
