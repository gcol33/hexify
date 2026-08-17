# Calculate resolution for target area

Uses the 'ISEA3H'/'ISEA4H'/'ISEA7H' cell count formula N = 10 \*
aperture^res + 2, which matches 'dggridR' resolution numbering exactly.

## Usage

``` r
calculate_resolution_for_area(
  target_area_km2,
  aperture = 3,
  radius_km = EARTH_RADIUS_KM
)
```

## Arguments

- target_area_km2:

  Target area in square kilometers

- aperture:

  Aperture (3, 4, or 7)

- radius_km:

  Radius of the body, in kilometers

## Value

Resolution level
