# Assign hexagonal DGGS cell IDs to geographic points

Takes a data.frame or sf object with geographic coordinates and returns
the same data with additional columns for hex cell ID and center
coordinates.

## Usage

``` r
hexify(
  data,
  lon = "lon",
  lat = "lat",
  area = NULL,
  diagonal = NULL,
  aperture = 3L,
  mixed_aperture_level = NULL,
  resround = "nearest"
)
```

## Arguments

- data:

  A data.frame or sf object containing coordinates

- lon:

  Column name for longitude (ignored if data is sf)

- lat:

  Column name for latitude (ignored if data is sf)

- area:

  Target cell area in km² (mutually exclusive with diagonal)

- diagonal:

  Target cell diagonal (long diagonal) in km (mutually exclusive with
  area)

- aperture:

  Grid aperture: 3, 4, 7, or "4/3" for mixed (default 3)

- mixed_aperture_level:

  For mixed aperture "4/3": number of aperture-4 levels before switching
  to aperture-3 (default NULL, auto-calculated as resolution/2)

- resround:

  How to round resolution: "nearest", "up", or "down"

## Value

The input data with additional columns:

- cell_id:

  Stable DGGS cell identifier (integer)

- cell_cen_lon:

  Longitude of cell center in degrees

- cell_cen_lat:

  Latitude of cell center in degrees

- cell_area:

  Actual cell area in km² (based on matched resolution)

- cell_diag:

  Actual cell long diagonal in km

## Details

For sf objects, coordinates are automatically extracted and transformed
to WGS84 (EPSG:4326) if needed. The geometry column is preserved.

Either `area` or `diagonal` must be provided, but not both:

- `area`: Target cell area in square kilometers

- `diagonal`: Long diagonal of hexagon in kilometers

The diagonal relates to area approximately as: area ≈ (3 \* sqrt(3) / 2)
\* (diagonal / 2)² ≈ 0.6495 \* diagonal²

Supported apertures:

- `aperture = 3`: ISEA3H (default, compatible with dggridR). Cell count
  = 10 \* 3^res + 2.

- `aperture = 4`: ISEA4H. Cell count = 10 \* 4^res + 2.

- `aperture = 7`: ISEA7H. Cell count = 10 \* 7^res + 2.

- `aperture = "4/3"`: ISEA43H mixed aperture. Uses aperture-4 for first
  `mixed_aperture_level` levels, then aperture-3. Cell count = 10 \*
  4^mixed_level \* 3^(res - mixed_level) + 2.

## See also

[`hexify_grid`](https://gcol33.github.io/hexify/reference/hexify_grid.md)
for grid specification,
[`hexify_to_sf`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)
for sf conversion

Other hexify main:
[`hexify_grid()`](https://gcol33.github.io/hexify/reference/hexify_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple data.frame
df <- data.frame(
  site = c("Vienna", "Paris", "Madrid"),
  lon = c(16.37, 2.35, -3.70),
  lat = c(48.21, 48.86, 40.42)
)
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# With sf object (any CRS)
library(sf)
pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
result_sf <- hexify(pts, area = 1000)

# Using diagonal instead of area
result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)

# Different apertures
result_ap4 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)
result_ap7 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)

# Mixed aperture (ISEA43H)
result_mixed <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "4/3")
} # }
```
