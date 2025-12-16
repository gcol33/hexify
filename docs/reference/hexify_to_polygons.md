# Generate polygons directly from hexify result

Convenience function that extracts resolution from a hexify result and
generates polygons. Resolution is auto-detected from the hex_area
column.

## Usage

``` r
hexify_to_polygons(data, aperture = 3L, return_sf = TRUE)
```

## Arguments

- data:

  Data frame returned by hexify() containing hex_id and hex_area

- aperture:

  Grid aperture (default 3)

- return_sf:

  Logical. If TRUE (default), returns sf object

## Value

sf object or data frame with polygon geometries (see hexify_cell_to_sf)

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
