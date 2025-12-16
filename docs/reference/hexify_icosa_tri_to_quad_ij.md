# Convert Icosa Triangle to Quad IJ coordinates

Converts icosahedral triangle coordinates directly to Quad IJ, combining
the transformation and quantization steps.

## Usage

``` r
hexify_icosa_tri_to_quad_ij(
  icosa_triangle_face,
  icosa_triangle_x,
  icosa_triangle_y,
  resolution,
  aperture = 3L
)
```

## Arguments

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

- resolution:

  Grid resolution level

- aperture:

  Grid aperture: 3, 4, or 7

## Value

List with components:

- quad:

  Quad number (0-11)

- i:

  Integer cell index along first axis

- j:

  Integer cell index along second axis

## Examples

``` r
if (FALSE) { # \dontrun{
# First get triangle coordinates from lon/lat
fwd <- hexify_forward(lon = 2.35, lat = 48.86)

# Then convert to quad IJ
quad_ij <- hexify_icosa_tri_to_quad_ij(
  icosa_triangle_face = fwd$face,
  icosa_triangle_x = fwd$tx,
  icosa_triangle_y = fwd$ty,
  resolution = 10,
  aperture = 3
)
print(quad_ij)
} # }
```
