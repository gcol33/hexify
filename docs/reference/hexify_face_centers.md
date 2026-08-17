# Get icosahedron face centers

Returns the center coordinates of all 20 icosahedral faces.

## Usage

``` r
hexify_face_centers()
```

## Value

Data frame with 20 rows and columns lon, lat (degrees)

## See also

Other projection:
[`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gillescolling.com/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gillescolling.com/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gillescolling.com/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gillescolling.com/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)

## Examples

``` r
centers <- hexify_face_centers()
plot(centers$lon, centers$lat)
```
