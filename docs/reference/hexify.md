# Assign hexagonal DGGS cell IDs to geographic points

Takes a data.frame or sf object with geographic coordinates and returns
a HexData object that stores the original data plus cell assignments.
The original data is preserved unchanged; cell IDs and centers are
stored in separate slots.

## Usage

``` r
hexify(
  data,
  grid = NULL,
  lon = "lon",
  lat = "lat",
  area_km2 = NULL,
  diagonal = NULL,
  resolution = NULL,
  aperture = 3,
  resround = "nearest",
  radius_km = EARTH_RADIUS_KM
)
```

## Arguments

- data:

  A data.frame or sf object containing coordinates

- grid:

  A HexGridInfo object from
  [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md).
  If provided, overrides area_km2, resolution, aperture and radius_km
  parameters.

- lon:

  Column name for longitude (ignored if data is sf)

- lat:

  Column name for latitude (ignored if data is sf)

- area_km2:

  Target cell area in km^2 (mutually exclusive with diagonal).

- diagonal:

  Target cell diagonal (long diagonal) in km

- resolution:

  Grid resolution (0-30). Alternative to area_km2.

- aperture:

  Grid aperture: 3, 4, 7, a mixed family such as "4/3" or "4/7", or one
  aperture per resolution level, e.g. `c(4, 4, 7, 3)` (default 3)

- resround:

  How to round resolution: "nearest", "up", or "down"

- radius_km:

  Radius of the body the grid covers, in kilometers, or a body name such
  as "mars" (default Earth). See
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md).

## Value

A HexData object containing:

- `data`: The original input data (unchanged)

- `grid`: The HexGridInfo specification

- `cell_id`: Numeric vector of cell IDs for each row

- `cell_center`: Matrix of cell center coordinates (lon, lat)

Use `as.data.frame(result)` to extract the original data. Use
`cells(result)` to get unique cell IDs. Use `result@cell_id` to get all
cell IDs. Use `result@cell_center` to get cell center coordinates.

## Details

For sf objects, coordinates are automatically extracted and transformed
to 'WGS84' (EPSG:4326) if needed. The geometry column is preserved.

Either `area_km2`, `diagonal`, or `resolution` must be provided unless a
`grid` object is supplied.

The HexData return type (default) stores the grid specification so
downstream functions like
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`hexify_cell_to_sf()`](https://gillescolling.com/hexify/reference/hexify_cell_to_sf.md),
etc. don't need grid parameters repeated.

## Grid Specification

You can create a grid specification once and reuse it:


    grid <- hex_grid(area_km2 = 1000)
    result1 <- hexify(df1, grid = grid)
    result2 <- hexify(df2, grid = grid)

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
grid specification,
[`HexData-class`](https://gillescolling.com/hexify/reference/HexData-class.md)
for return object details,
[`as_sf`](https://gillescolling.com/hexify/reference/as_sf.md) for
converting to sf

Other hexify main:
[`hexify_grid()`](https://gillescolling.com/hexify/reference/hexify_grid.md)

## Examples

``` r
# Simple data.frame
df <- data.frame(
  site = c("Vienna", "Paris", "Madrid"),
  lon = c(16.37, 2.35, -3.70),
  lat = c(48.21, 48.86, 40.42)
)

# New recommended workflow: use grid object
grid <- hex_grid(area_km2 = 1000)
result <- hexify(df, grid = grid, lon = "lon", lat = "lat")
print(result)  # Shows grid info
plot(result)   # Plot with default styling

# Direct area specification (grid created internally)
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Extract plain data.frame
df_result <- as.data.frame(result)

# With sf object (any CRS)
library(sf)
pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
result_sf <- hexify(pts, area_km2 = 1000)

# Different apertures
result_ap4 <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 4)

# Mixed aperture (ISEA43H)
result_mixed <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = "4/3")
```
