# Forward Snyder projection

Projects geographic coordinates onto the icosahedron, returning face
index and planar coordinates (tx, ty).

## Usage

``` r
hexify_forward(lon, lat)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

## Value

Named numeric vector: c(face, tx, ty)

## Details

tx and ty are normalized coordinates within the triangular face,
typically in range 0, 1.

## Examples

``` r
result <- hexify_forward(16.37, 48.21)
# result["face"], result["icosa_triangle_x"], result["icosa_triangle_y"]
```
