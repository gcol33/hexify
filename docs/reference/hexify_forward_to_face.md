# Forward projection to specific face

Projects to a known face (skips face detection).

## Usage

``` r
hexify_forward_to_face(face, lon, lat)
```

## Arguments

- face:

  Face index (0-19)

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

## Value

Named numeric vector: c(icosa_triangle_x, icosa_triangle_y)

## See also

Other projection:
[`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md),
[`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gcol33.github.io/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)
