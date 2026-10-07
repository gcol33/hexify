# Forward face projection

Projects geographic coordinates onto a face of the solid, returning face
index and planar coordinates (tx, ty).

## Usage

``` r
hexify_forward(
  lon,
  lat,
  projection = c("isea", "fuller"),
  polyhedron = c("icosahedron", "octahedron", "tetrahedron")
)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

- projection:

  Face projection: `"isea"` (Snyder's equal-area projection) or
  `"fuller"` (Fuller's projection, defined on the icosahedron only)

- polyhedron:

  The solid: "icosahedron" (default), "octahedron" or "tetrahedron", in
  its default orientation. Snyder (1992) gives his equal-area projection
  for each; its angular distortion grows with the faces' size.

## Value

Named numeric vector: c(face, tx, ty)

## Details

tx and ty are normalized coordinates within the triangular face,
typically in range \[0, 1\].

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)

## Examples

``` r
result <- hexify_forward(16.37, 48.21)
# result["face"], result["icosa_triangle_x"], result["icosa_triangle_y"]
hexify_forward(16.37, 48.21, polyhedron = "tetrahedron")
```
