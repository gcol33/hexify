# Find closest resolution for target cell spacing

Finds the grid resolution that produces cells with spacing (distance
between centers) closest to the target spacing.

## Usage

``` r
dg_closest_res_to_spacing(
  dggs,
  spacing,
  round = "nearest",
  metric = TRUE,
  show_info = FALSE
)
```

## Arguments

- dggs:

  Grid specification (aperture and topology must be set)

- spacing:

  Target cell spacing in km (if metric=TRUE)

- round:

  Rounding method ("nearest", "up", "down")

- metric:

  Whether spacing is in metric units

- show_info:

  Print information about chosen resolution

## Value

Resolution level (integer)

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`dg_closest_res_to_cls()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_cls.md),
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)
