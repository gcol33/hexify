# Compare grid resolutions

Generates a table comparing different resolution levels for a given grid
configuration. Useful for choosing appropriate resolution.

## Usage

``` r
hexify_compare_resolutions(aperture = 3, res_range = 0:15)
```

## Arguments

- aperture:

  Grid aperture (3, 4, or 7)

- res_range:

  Range of resolutions to compare (e.g., 1:10)

## Value

Data frame with columns: resolution, n_cells, cell_area_km2,
cell_spacing_km, cls_km

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`dg_closest_res_to_cls()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_cls.md),
[`dg_closest_res_to_spacing()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_spacing.md),
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare resolutions 0-10 for aperture 3
comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:10)
print(comparison)

# Find resolution with cells ~1000 km²
subset(comparison, cell_area_km2 > 900 & cell_area_km2 < 1100)
} # }
```
