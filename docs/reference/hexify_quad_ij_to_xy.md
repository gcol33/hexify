# Convert Quad IJ to Quad XY (continuous coordinates)

Converts discrete cell indices to continuous quad coordinates. Useful
for computing cell centers or understanding the cell geometry.

## Usage

``` r
hexify_quad_ij_to_xy(quad, i, j, resolution, aperture = 3L)
```

## Arguments

- quad:

  Quad number (0-11)

- i:

  Cell index along first axis

- j:

  Cell index along second axis

- resolution:

  Grid resolution level

- aperture:

  Grid aperture: 3, 4, or 7

## Value

List with components:

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## Examples

``` r
if (FALSE) { # \dontrun{
# Get continuous quad coordinates for a cell
xy <- hexify_quad_ij_to_xy(quad = 1, i = 100, j = 50,
                           resolution = 10, aperture = 3)
print(xy)
} # }
```
