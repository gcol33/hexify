# Generate polygons directly from hexify result

Convenience function that extracts resolution from a hexify result and
generates polygons. Resolution is auto-detected from the cell_area
column.

## Usage

``` r
hexify_to_polygons(data, aperture = 3L, return_sf = TRUE)
```

## Arguments

- data:

  Data frame returned by hexify() containing cell_id and cell_area

- aperture:

  Grid aperture (default 3)

- return_sf:

  Logical. If TRUE (default), returns sf object

## Value

sf object or data frame with polygon geometries (see hexify_cell_to_sf)

## See also

[`hexify_cell_to_sf`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md)
for low-level polygon generation,
[`hexify_to_sf`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)
for full hexify result conversion

Other sf conversion:
[`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md),
[`hexify_cell_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md),
[`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md),
[`hexify_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)
library(sf)

# Simple workflow: hexify then get polygons
df <- data.frame(lon = runif(100, -10, 10), lat = runif(100, 40, 50))
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# Get polygons - resolution auto-detected
polys <- hexify_to_polygons(result)
plot(st_geometry(polys))
} # }
```
