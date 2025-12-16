# Build an sf POLYGON from six (lon, lat) corner pairs

Build an sf POLYGON from six (lon, lat) corner pairs

## Usage

``` r
hex_corners_to_sf(lon, lat, crs = 4326)
```

## Arguments

- lon:

  numeric vector of length 6 (longitude)

- lat:

  numeric vector of length 6 (latitude)

- crs:

  integer CRS (default 4326)

## Value

sf object with one POLYGON geometry
