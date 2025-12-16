# Convert Quad IJ to Icosa Triangle coordinates

Converts Quad IJ coordinates to icosahedral triangle coordinates. This
is useful for understanding where a cell is located on the icosahedral
projection.

## Usage

``` r
hexify_quad_ij_to_icosa_tri(quad, i, j, resolution, aperture = 3L)
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

Data frame with columns:

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Details

Equivalent to dggridR's dgQ2DI_to_PROJTRI().

## See also

[`hexify_icosa_tri_to_quad_ij`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_quad_ij.md)
for the inverse,
[`hexify_cell_to_icosa_tri`](https://gcol33.github.io/hexify/reference/hexify_cell_to_icosa_tri.md)
for conversion from cell ID

## Examples

``` r
if (FALSE) { # \dontrun{
# Get triangle coordinates for a Quad IJ position
result <- hexify_quad_ij_to_icosa_tri(quad = 1, i = 100, j = 50,
                                       resolution = 10, aperture = 3)
print(result)
} # }
```
