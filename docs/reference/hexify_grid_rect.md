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
  resround = "nearest",
  radius_km = EARTH_RADIUS_KM
)
```

## Arguments

- minlon, maxlon:

  Longitude bounds

- minlat, maxlat:

  Latitude bounds

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

sf object with hexagon polygons covering the specified region

## See also

[`grid_rect`](https://gillescolling.com/hexify/reference/grid_rect.md)
for the recommended S4 interface,
[`hexify_grid_global`](https://gillescolling.com/hexify/reference/hexify_grid_global.md)
for global grids

Other sf conversion:
[`hex_corners_to_sf()`](https://gillescolling.com/hexify/reference/hex_corners_to_sf.md),
[`hexify_cell_to_sf()`](https://gillescolling.com/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_global()`](https://gillescolling.com/hexify/reference/hexify_grid_global.md)

## Examples

``` r
library(hexify)
library(sf)

grid <- hexify_grid_rect(
  minlon = -10, maxlon = 20,
  minlat = 35, maxlat = 60,
  area = 5000
)
plot(st_geometry(grid), border = "gray")
```
