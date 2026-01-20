# Convert HexData to sf Object

Converts a HexData object to an sf spatial features object. Can create
either point geometries (cell centers) or polygon geometries (cell
boundaries).

## Usage

``` r
as_sf(x, geometry = c("point", "polygon"), ...)
```

## Arguments

- x:

  A HexData object

- geometry:

  Type of geometry: "point" (default) or "polygon"

- ...:

  Additional arguments (ignored)

## Value

An sf object

## Details

For point geometry, cell centers (cell_cen_lon, cell_cen_lat) are used.
For polygon geometry, cell boundaries are computed using the grid
specification.

## Examples

``` r
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Get sf points
sf_pts <- as_sf(result)

# Get sf polygons
sf_poly <- as_sf(result, geometry = "polygon")
```
