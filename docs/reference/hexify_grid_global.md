# Generate a global grid of hexagon polygons

Creates hexagon polygons covering the entire Earth.

## Usage

``` r
hexify_grid_global(area, aperture = 3L, resround = "nearest")
```

## Arguments

- area:

  Target cell area in km²

- aperture:

  Grid aperture: 3, 4, or 7

- resround:

  Resolution rounding: "nearest", "up", or "down"

## Value

sf object with hexagon polygons covering the globe

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)
library(sf)

# Coarse global grid (~100,000 km² cells)
global_grid <- hexify_grid_global(area = 100000)
plot(st_geometry(global_grid), border = "gray")
} # }
```
