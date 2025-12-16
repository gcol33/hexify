# Convert Cell ID to Quad IJ coordinates

Converts DGGRID-compatible cell IDs to Quad IJ coordinates. This is the
inverse of hexify_quad_ij_to_cell().

## Usage

``` r
hexify_cell_to_quad_ij(cell_id, resolution, aperture = 3L)
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

- i:

  Integer cell index along first axis

- j:

  Integer cell index along second axis

## Details

Compatible with dggridR's dgSEQNUM_to_Q2DI().

## See also

[`hexify_quad_ij_to_cell`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md)
for the forward operation,
[`hexify_cell_to_icosa_tri`](https://gcol33.github.io/hexify/reference/hexify_cell_to_icosa_tri.md)
for conversion to triangle coords

## Examples

``` r
if (FALSE) { # \dontrun{
# Get Quad IJ coordinates for a cell
result <- hexify_cell_to_quad_ij(cell_id = 1000, resolution = 10, aperture = 3)
print(result)

# Round-trip test
cell_id <- hexify_quad_ij_to_cell(result$quad, result$i, result$j,
                                   resolution = 10, aperture = 3)
# Should equal original cell_id
} # }
```
