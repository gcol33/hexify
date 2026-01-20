# Compare grid resolutions

Generates a table comparing different resolution levels for a given grid
configuration. Useful for choosing appropriate resolution.

## Usage

``` r
hexify_compare_resolutions(aperture = 3, res_range = 0:15, print = FALSE)
```

## Arguments

- aperture:

  Grid aperture (3, 4, or 7)

- res_range:

  Range of resolutions to compare (e.g., 1:10)

- print:

  If TRUE, prints a formatted table to console. If FALSE (default),
  returns a data frame.

## Value

If print=FALSE: data frame with columns resolution, n_cells,
cell_area_km2, cell_spacing_km, cls_km. If print=TRUE: invisibly returns
the data frame after printing.

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)

## Examples

``` r
# Get data frame of resolutions 0-10 for aperture 3
comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:10)
print(comparison)

# Print formatted table directly
hexify_compare_resolutions(aperture = 3, res_range = 0:10, print = TRUE)

# Find resolution with cells ~1000 km^2
subset(comparison, cell_area_km2 > 900 & cell_area_km2 < 1100)
```
