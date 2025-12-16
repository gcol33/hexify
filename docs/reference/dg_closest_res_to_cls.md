# Find closest resolution for target CLS

Finds the grid resolution that produces cells with characteristic length
scale (CLS) closest to the target CLS.

## Usage

``` r
dg_closest_res_to_cls(
  dggs,
  cls,
  round = "nearest",
  metric = TRUE,
  show_info = FALSE
)
```

## Arguments

- dggs:

  Grid specification (aperture and topology must be set)

- cls:

  Target CLS in km (if metric=TRUE)

- round:

  Rounding method ("nearest", "up", "down")

- metric:

  Whether CLS is in metric units

- show_info:

  Print information about chosen resolution

## Value

Resolution level (integer)
