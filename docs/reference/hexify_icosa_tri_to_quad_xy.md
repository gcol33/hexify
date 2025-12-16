# Convert Icosa Triangle to Quad XY coordinates

Converts icosahedral triangle coordinates (from Snyder projection) to
quad XY coordinates. This is an intermediate step in the pipeline.

## Usage

``` r
hexify_icosa_tri_to_quad_xy(
  icosa_triangle_face,
  icosa_triangle_x,
  icosa_triangle_y
)
```

## Arguments

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Value

List with components:

- quad:

  Quad number (0-11)

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## Examples

``` r
if (FALSE) { # \dontrun{
# First get triangle coordinates from lon/lat
fwd <- hexify_forward(lon = 2.35, lat = 48.86)

# Then convert to quad XY
quad_xy <- hexify_icosa_tri_to_quad_xy(
  icosa_triangle_face = fwd$face,
  icosa_triangle_x = fwd$tx,
  icosa_triangle_y = fwd$ty
)
print(quad_xy)
} # }
```
