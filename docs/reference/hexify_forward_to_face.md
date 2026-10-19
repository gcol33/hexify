# Forward Projection to Specific Face

Projects to a known face (skips face detection).

## Usage

``` r
hexify_forward_to_face(
  face,
  lon,
  lat,
  projection = c("isea", "fuller", "ivea", "ak", "akw"),
  polyhedron = c("icosahedron", "octahedron", "tetrahedron")
)
```

## Arguments

- face:

  Face index (0-19 on the icosahedron)

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

- projection:

  Face projection: `"isea"` (Snyder's equal-area projection), `"ivea"`
  (van Leeuwen and Strebe's vertex-oriented equal-area projection),
  `"fuller"` (Fuller's projection, defined on the icosahedron only),
  `"ak"` (Kaseorg's octahedral projection) or `"akw"` (Kaseorg's with
  Hex9's trained warp), both defined on the octahedron only

- polyhedron:

  The solid: "icosahedron" (default), "octahedron" or "tetrahedron", in
  its default orientation. Snyder (1992) gives his equal-area projection
  for each; its angular distortion grows with the faces' size.

## Value

Named numeric vector: c(icosa_triangle_x, icosa_triangle_y)

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)
