# Determine which face contains a point

Returns the icosahedral face index (0-19) containing the given
coordinates.

## Usage

``` r
hexify_which_face(lon, lat)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

## Value

Integer face index (0-19)

## See also

Other projection:
[`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gcol33.github.io/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md)

## Examples

``` r
face <- hexify_which_face(16.37, 48.21)
```
