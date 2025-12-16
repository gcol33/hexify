# Convert longitude/latitude to index string

Main entry point for geocoding points to grid cells.

## Usage

``` r
hexify_lonlat_to_index(
  lon,
  lat,
  resolution,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

- resolution:

  Resolution level

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding: "auto" (default), "z3", "z7", or "zorder"

## Value

Index string

## Examples

``` r
idx <- hexify_lonlat_to_index(16.37, 48.21, resolution = 5, aperture = 3)
```
