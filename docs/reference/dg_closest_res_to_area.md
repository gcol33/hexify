# Find closest resolution for target cell area

Finds the grid resolution that produces cells closest to the target
area. This is primarily used internally by
[`hexify_grid`](https://gcol33.github.io/hexify/reference/hexify_grid.md)
and [`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md).
Most users should use those functions directly.

## Usage

``` r
dg_closest_res_to_area(
  dggs,
  area,
  round = "nearest",
  metric = TRUE,
  show_info = FALSE
)
```

## Arguments

- dggs:

  Grid specification (aperture and topology must be set)

- area:

  Target cell area in km^2 (if metric=TRUE)

- round:

  Rounding method ("nearest", "up", "down")

- metric:

  Whether area is in metric units

- show_info:

  Print information about chosen resolution

## Value

Resolution level (integer)

## See also

Other grid statistics:
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a temporary grid to get aperture settings
temp_grid <- list(aperture = 3, topology = "HEXAGON")
class(temp_grid) <- "hexify_grid"

# Find resolution for 1000 km^2 cells
res <- dg_closest_res_to_area(temp_grid, area = 1000, 
                               metric = TRUE, show_info = TRUE)
print(res)
} # }
```
