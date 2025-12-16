# Convert Quad IJ coordinates to cell ID

Converts Quad IJ coordinates to a global cell identifier. This is the
final step in the coordinate pipeline.

## Usage

``` r
hexify_quad_ij_to_cell(quad, i, j, resolution, aperture = 3L)
```

## Arguments

- quad:

  Quad number (0-11), integer or vector

- i:

  Cell index along first axis, integer or vector

- j:

  Cell index along second axis, integer or vector

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Numeric vector of cell IDs

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Quad IJ to cell ID
cell_id <- hexify_quad_ij_to_cell(quad = 1, i = 100, j = 50,
                                   resolution = 10, aperture = 3)
print(cell_id)
} # }
```
