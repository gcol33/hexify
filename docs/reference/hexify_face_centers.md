# Get face centers of a solid

Returns the center coordinates of every face of the solid: 20 on the
icosahedron, 8 on the octahedron, 4 on the tetrahedron.

## Usage

``` r
hexify_face_centers(polyhedron = c("icosahedron", "octahedron", "tetrahedron"))
```

## Arguments

- polyhedron:

  The solid: "icosahedron" (default), "octahedron" or "tetrahedron", in
  its default orientation

## Value

Data frame with one row per face and columns lon, lat (radians)

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)

## Examples

``` r
centers <- hexify_face_centers()
plot(centers$lon, centers$lat)
nrow(hexify_face_centers("octahedron"))
```
