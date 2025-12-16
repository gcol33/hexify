# Convert cell ID to longitude/latitude

Converts cell identifiers back to cell center coordinates. This is the
inverse of
[`hexify_lonlat_to_cell`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md).

## Usage

``` r
hexify_cell_to_lonlat(cell_id, resolution, aperture)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs (1-based)

- resolution:

  Grid resolution (integer \>= 0)

- aperture:

  Grid aperture (3, 4, or 7)

## Value

Data frame with lon_deg and lat_deg columns

## See also

[`hexify_lonlat_to_cell`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md)
for the forward operation,
[`hexify_grid_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md)
for the grid-based wrapper

## Examples

``` r
coords <- hexify_cell_to_lonlat(1702, resolution = 5, aperture = 3)
```
