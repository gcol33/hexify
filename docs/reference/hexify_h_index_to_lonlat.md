# Convert hierarchical index strings to longitude/latitude centers

Converts hierarchical cell index strings back to geographic coordinates,
returning the center point of each cell. This is the inverse operation
of hexify_lonlat_to_h_index().

## Usage

``` r
hexify_h_index_to_lonlat(grid, h_index)
```

## Arguments

- grid:

  Grid specification from hexify_grid()

- h_index:

  Hierarchical index strings (character vector)

## Value

Data frame with columns:

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

## Details

Most users should use
[`hexify_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md)
or
[`hexify_grid_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md)
which work with DGGRID-compatible integer cell IDs.
