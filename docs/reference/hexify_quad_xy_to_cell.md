# Convert Quad XY coordinates to Cell ID

Converts Quad XY coordinates (continuous quad space) to
DGGRID-compatible cell IDs. The coordinates are quantized to the nearest
cell.

## Usage

``` r
hexify_quad_xy_to_cell(quad, quad_x, quad_y, resolution, aperture = 3L)
```

## Arguments

- quad:

  Quad number (0-11), integer or vector

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Numeric vector of cell IDs

## Details

Compatible with dggridR's dgQ2DD_to_SEQNUM().

## See also

[`hexify_cell_to_quad_xy`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_xy.md)
for the inverse operation,
[`hexify_quad_ij_to_cell`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md)
for integer grid coordinates

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Quad XY to cell ID
cell_id <- hexify_quad_xy_to_cell(quad = 1, quad_x = 0.5, quad_y = 0.3,
                                   resolution = 10, aperture = 3)
print(cell_id)
} # }
```
