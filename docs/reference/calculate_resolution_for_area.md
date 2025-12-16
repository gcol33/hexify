# Calculate resolution for target area

Uses the ISEA3H cell count formula: N = 10 \* aperture^res + 2 This
matches dggridR's resolution numbering exactly.

## Usage

``` r
calculate_resolution_for_area(target_area_km2, aperture = 3)
```

## Arguments

- target_area_km2:

  Target area in square kilometers

- aperture:

  Aperture (3, 4, or 7)

## Value

Resolution level
