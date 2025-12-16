# Convert Cell ID to Quad XY coordinates

Converts DGGRID-compatible cell IDs to Quad XY coordinates (continuous
quad space). This is the cell center in quad coordinates.

## Usage

``` r
hexify_cell_to_quad_xy(cell_id, resolution, aperture = 3L)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs (1-based)

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Data frame with columns:

- quad:

  Quad number (0-11)

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## Details

Compatible with dggridR's dgSEQNUM_to_Q2DD().

## See also

[`hexify_quad_xy_to_cell`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md)
for the inverse operation,
[`hexify_cell_to_quad_ij`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_ij.md)
for integer grid coordinates

## Examples

``` r
if (FALSE) { # \dontrun{
# Get Quad XY coordinates for a cell
result <- hexify_cell_to_quad_xy(cell_id = 1000, resolution = 10, aperture = 3)
print(result)

# Round-trip test
cell_id <- hexify_quad_xy_to_cell(result$quad, result$quad_x, result$quad_y,
                                   resolution = 10, aperture = 3)
# Should equal original cell_id
} # }
```
