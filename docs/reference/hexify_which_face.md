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
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gillescolling.com/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gillescolling.com/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gillescolling.com/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gillescolling.com/hexify/reference/hexify_set_verbose.md)

## Examples

``` r
face <- hexify_which_face(16.37, 48.21)
```
