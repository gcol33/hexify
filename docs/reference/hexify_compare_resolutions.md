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
