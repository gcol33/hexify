# Convert longitude/latitude to cell ID using a grid object

Grid-based wrapper for
[`hexify_lonlat_to_cell`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md).
Converts geographic coordinates to DGGRID-compatible cell IDs using the
resolution and aperture from a grid object.

## Usage

``` r
hexify_grid_to_cell(grid, lon, lat)
```

## Arguments

- grid:

  Grid specification from hexify_grid()

- lon:

  Numeric vector of longitudes in degrees

- lat:

  Numeric vector of latitudes in degrees

## Value

Numeric vector of cell IDs (1-based)

## See also

[`hexify_lonlat_to_cell`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md)
for the direct-params version,
[`hexify_grid_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md)
for the inverse operation

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hexify_grid(area = 1000, aperture = 3)
cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))
} # }
```
