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
