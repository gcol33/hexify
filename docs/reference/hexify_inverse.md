# Inverse face projection

Converts face plane coordinates back to geographic coordinates.

## Usage

``` r
hexify_inverse(
  x,
  y,
  face,
  projection = c("isea", "fuller"),
  polyhedron = c("icosahedron", "octahedron", "tetrahedron")
)
```

## Arguments

- x:

  X coordinate on face plane

- y:

  Y coordinate on face plane

- face:

  Face index (0-19 on the icosahedron)

- projection:

  Face projection: `"isea"` (Snyder's equal-area projection) or
  `"fuller"` (Fuller's projection, defined on the icosahedron only)

- polyhedron:

  The solid: "icosahedron" (default), "octahedron" or "tetrahedron", in
  its default orientation. Snyder (1992) gives his equal-area projection
  for each; its angular distortion grows with the faces' size.

## Value

Named numeric vector: c(lon_deg, lat_deg)

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)

## Examples

``` r
coords <- hexify_inverse(0.5, 0.3, face = 2)
```
