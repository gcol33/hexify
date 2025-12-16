# Convert longitude/latitude to Quad IJ coordinates

Converts geographic coordinates to the intermediate Quad IJ
representation used internally by ISEA DGGS. Returns the quad number
(0-11) and integer cell indices (i, j) within that quad.

## Usage

``` r
hexify_lonlat_to_quad_ij(lon, lat, resolution, aperture = 3L)
```

## Arguments

- lon:

  Longitude in degrees (-180 to 180)

- lat:

  Latitude in degrees (-90 to 90)

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

List with components:

- quad:

  Quad number (0-11)

- i:

  Integer cell index along first axis

- j:

  Integer cell index along second axis

- icosa_triangle_face:

  Source icosahedral face (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Details

The 20 icosahedral triangle faces are grouped into 12 quads:

- Quad 0: North polar region

- Quads 1-5: Upper hemisphere rhombi

- Quads 6-10: Lower hemisphere rhombi

- Quad 11: South polar region

## Examples

``` r
if (FALSE) { # \dontrun{
# Get Quad IJ coordinates for Paris
result <- hexify_lonlat_to_quad_ij(lon = 2.35, lat = 48.86,
                                    resolution = 10, aperture = 3)
print(result)
} # }
```
