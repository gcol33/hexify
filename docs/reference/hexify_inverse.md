# Inverse Snyder projection

Converts face plane coordinates back to geographic coordinates.

## Usage

``` r
hexify_inverse(x, y, face, tol = NULL, max_iters = NULL)
```

## Arguments

- x:

  X coordinate on face plane

- y:

  Y coordinate on face plane

- face:

  Face index (0-19)

- tol:

  Convergence tolerance (NULL for default)

- max_iters:

  Maximum iterations (NULL for default)

## Value

Named numeric vector: c(lon_deg, lat_deg)

## See also

Other projection:
[`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md),
[`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gcol33.github.io/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)

## Examples

``` r
coords <- hexify_inverse(0.5, 0.3, face = 2)
```
