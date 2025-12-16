# Round-trip accuracy test

Tests the accuracy of the coordinate conversion functions by converting
coordinates to cells and back, measuring the distance between original
and reconstructed coordinates.

## Usage

``` r
hexify_roundtrip_test(grid, lon, lat, units = "km")
```

## Arguments

- grid:

  Grid specification

- lon:

  Longitude to test

- lat:

  Latitude to test

- units:

  Distance units ("km" or "degrees")

## Value

List with:

- original:

  Original coordinates

- cell:

  Cell index

- reconstructed:

  Reconstructed coordinates

- error:

  Distance between original and reconstructed

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hexify_grid(area = 1000, aperture = 3)
accuracy <- hexify_roundtrip_test(grid, lon = 0, lat = 45)
print(accuracy)
} # }
```
