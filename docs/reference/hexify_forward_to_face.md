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
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_get_precision()`](https://gillescolling.com/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gillescolling.com/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gillescolling.com/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gillescolling.com/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)
