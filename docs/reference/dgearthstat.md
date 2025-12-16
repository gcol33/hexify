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

  Total Earth surface area in km²

- n_cells:

  Total number of cells at this resolution

- cell_area_km2:

  Average cell area in km²

- cell_spacing_km:

  Average distance between cell centers in km

- resolution:

  Resolution level

- aperture:

  Grid aperture

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hexify_grid(area = 1000, aperture = 3)
stats <- dgearthstat(grid)

print(sprintf("Resolution %d has %.0f cells",
              stats$resolution, stats$n_cells))
print(sprintf("Average cell area: %.2f km²",
              stats$cell_area_km2))
print(sprintf("Average cell spacing: %.2f km",
              stats$cell_spacing_km))
} # }
```
