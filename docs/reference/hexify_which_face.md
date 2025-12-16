# Determine which face contains a point

Returns the icosahedral face index (0-19) containing the given
coordinates.

## Usage

``` r
hexify_which_face(lon, lat)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

## Value

Integer face index (0-19)

## Examples

``` r
face <- hexify_which_face(16.37, 48.21)
```
