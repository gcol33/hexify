# Find closest resolution for target cell area

Finds the grid resolution that produces cells closest to the target
area. This is a helper function for grid construction.

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

  Target cell area in km² (if metric=TRUE)

- round:

  Rounding method ("nearest", "up", "down")

- metric:

  Whether area is in metric units

- show_info:

  Print information about chosen resolution

## Value

Resolution level (integer)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a temporary grid to get aperture settings
temp_grid <- list(aperture = 3, topology = "HEXAGON")
class(temp_grid) <- "hexify_grid"

# Find resolution for 1000 km² cells
res <- dg_closest_res_to_area(temp_grid, area = 1000, 
                               metric = TRUE, show_info = TRUE)
print(res)
} # }
```
