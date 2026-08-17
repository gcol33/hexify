# Generate a global grid of hexagon polygons

Creates hexagon polygons covering a whole body.

## Usage

``` r
hexify_grid_global(
  area,
  aperture = 3L,
  resround = "nearest",
  radius_km = EARTH_RADIUS_KM
)
```

## Arguments

- area:

  Target cell area in km^2

- aperture:

  Grid aperture: 3, 4, or 7

- resround:

  Resolution rounding: "nearest", "up", or "down"

- radius_km:

  Radius of the body, in kilometers, or a body name such as "mars"
  (default Earth). See
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md).

## Value

sf object with hexagon polygons covering the globe

## See also

[`grid_global`](https://gillescolling.com/hexify/reference/grid_global.md)
for the recommended S4 interface,
[`hexify_grid_rect`](https://gillescolling.com/hexify/reference/hexify_grid_rect.md)
for regional grids

Other sf conversion:
[`hex_corners_to_sf()`](https://gillescolling.com/hexify/reference/hex_corners_to_sf.md),
[`hexify_cell_to_sf()`](https://gillescolling.com/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_rect()`](https://gillescolling.com/hexify/reference/hexify_grid_rect.md)

## Examples

``` r
library(hexify)
library(sf)

# Coarse global grid (~100,000 km^2 cells)
global_grid <- hexify_grid_global(area = 100000)
plot(st_geometry(global_grid), border = "gray")
```
