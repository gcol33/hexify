# Get icosahedron face centers

Returns the center coordinates of all 20 icosahedral faces.

## Usage

``` r
hexify_face_centers()
```

## Value

Data frame with 20 rows and columns lon, lat (degrees)

## Examples

``` r
centers <- hexify_face_centers()
plot(centers$lon, centers$lat)
```
