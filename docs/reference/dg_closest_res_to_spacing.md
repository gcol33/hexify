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
