# Determine which face contains a point

Returns the index of the face of the solid containing the given
coordinates (0-19 on the icosahedron).

## Usage

``` r
hexify_which_face(
  lon,
  lat,
  polyhedron = c("icosahedron", "octahedron", "tetrahedron")
)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

- polyhedron:

  The solid: "icosahedron" (default), "octahedron" or "tetrahedron", in
  its default orientation

## Value

Integer face index

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md)

## Examples

``` r
face <- hexify_which_face(16.37, 48.21)
hexify_which_face(16.37, 48.21, polyhedron = "tetrahedron")
```
