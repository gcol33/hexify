# Get cell info from cell ID

Converts cell ID to cell components (quad, i, j).

## Usage

``` r
hexify_cell_id_to_quad_ij(cell_id, resolution, aperture)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs (1-based)

- resolution:

  Grid resolution (integer \>= 0)

- aperture:

  Grid aperture (3, 4, or 7)

## Value

Data frame with quad, i, j columns

## Examples

``` r
info <- hexify_cell_to_quad_ij(1702, resolution = 5, aperture = 3)
```
