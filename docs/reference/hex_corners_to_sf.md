# Build an sf POLYGON from six (lon, lat) corner pairs

Build an sf POLYGON from six (lon, lat) corner pairs

## Usage

``` r
hex_corners_to_sf(lon, lat, crs = 4326)
```

## Arguments

- lon:

  numeric vector of length 6 (longitude)

- lat:

  numeric vector of length 6 (latitude)

- crs:

  integer CRS (default 4326)

## Value

sf object with one POLYGON geometry

## See also

Other sf conversion:
[`hexify_cell_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md),
[`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md),
[`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md),
[`hexify_to_polygons()`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md),
[`hexify_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)
