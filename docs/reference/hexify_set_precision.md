# Set inverse projection precision

Controls the accuracy/speed tradeoff for inverse Snyder projection.

## Usage

``` r
hexify_set_precision(
  mode = c("fast", "default", "high", "ultra"),
  tol = NULL,
  max_iters = NULL
)
```

## Arguments

- mode:

  Preset mode: "fast", "default", "high", or "ultra"

- tol:

  Custom tolerance (overrides mode if provided)

- max_iters:

  Custom max iterations (overrides mode if provided)

## Value

Invisible NULL

## See also

Other projection:
[`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md),
[`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md),
[`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)

## Examples

``` r
hexify_set_precision("high")
hexify_set_precision(tol = 1e-12, max_iters = 100)
```
