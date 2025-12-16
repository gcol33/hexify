# Convert cell ID to longitude/latitude using a grid object

Grid-based wrapper for
[`hexify_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md).
Converts DGGRID-compatible cell IDs back to cell center coordinates
using the resolution and aperture from a grid object.

## Usage

``` r
hexify_grid_cell_to_lonlat(grid, cell_id)
```

## Arguments

- grid:

  Grid specification from hexify_grid()

- cell_id:

  Numeric vector of cell IDs (1-based)

## Value

Data frame with lon_deg and lat_deg columns

## See also

[`hexify_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md)
for the direct-params version,
[`hexify_grid_to_cell`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md)
for the forward operation

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hexify_grid(area = 1000, aperture = 3)
cell_ids <- hexify_grid_to_cell(grid, lon = 5, lat = 45)
coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)
} # }
```
