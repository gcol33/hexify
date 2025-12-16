# Convert Quad XY to Icosa Triangle coordinates

Inverse transformation from quad coordinates back to icosahedral
triangle coordinates. Useful for projecting cell centers back to
lon/lat.

## Usage

``` r
hexify_quad_xy_to_icosa_tri(quad, quad_x, quad_y)
```

## Arguments

- quad:

  Quad number (0-11)

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## Value

List with components:

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert quad XY back to triangle coordinates
tri <- hexify_quad_xy_to_icosa_tri(quad = 1, quad_x = 0.5, quad_y = 0.3)
print(tri)
} # }
```
