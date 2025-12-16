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

## Examples

``` r
hexify_set_precision("high")
hexify_set_precision(tol = 1e-12, max_iters = 100)
```
