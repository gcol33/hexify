# Get grid statistics for Earth coverage

Calculates statistics about the hexagonal grid at the current
resolution, including total number of cells, cell area, and cell
spacing.

## Usage

``` r
dgearthstat(dggs)
```

## Arguments

- dggs:

  Grid specification from hexify_grid()

## Value

List with components:

- area_km:

  Total Earth surface area in km^2

- n_cells:

  Total number of cells at this resolution

- cell_area_km2:

  Average cell area in km^2

- cell_spacing_km:

  Average distance between cell centers in km

- resolution:

  Resolution level

- aperture:

  Grid aperture

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hexify_grid(area = 1000, aperture = 3)
stats <- dgearthstat(grid)

print(sprintf("Resolution %d has %.0f cells",
              stats$resolution, stats$n_cells))
print(sprintf("Average cell area: %.2f km^2",
              stats$cell_area_km2))
print(sprintf("Average cell spacing: %.2f km",
              stats$cell_spacing_km))
} # }
```
