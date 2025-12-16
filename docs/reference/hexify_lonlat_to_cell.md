# Convert longitude/latitude to cell ID

Converts geographic coordinates to DGGRID-compatible cell identifiers.
This is the primary function for geocoding points to grid cells.

## Usage

``` r
hexify_lonlat_to_cell(lon, lat, resolution, aperture)
```

## Arguments

- lon:

  Numeric vector of longitudes in degrees

- lat:

  Numeric vector of latitudes in degrees

- resolution:

  Grid resolution (integer \>= 0)

- aperture:

  Grid aperture (3, 4, or 7)

## Value

Numeric vector of cell IDs (1-based)

## Details

Returns DGGRID-compatible cell identifiers. The cell ID uniquely
identifies each hexagonal cell in the global grid.

For a grid-object interface, use
[`hexify_grid_to_cell`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md).

## See also

[`hexify_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md)
for the inverse operation,
[`hexify_grid_to_cell`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md)
for the grid-based wrapper

## Examples

``` r
cell_id <- hexify_lonlat_to_cell(0, 45, resolution = 5, aperture = 3)
```
