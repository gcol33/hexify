# Assign hexagonal DGGS cell IDs to geographic points

Takes a data.frame or sf object with geographic coordinates and returns
the data with additional columns for hex cell ID and center coordinates.
By default returns a HexData object that stores the grid specification
for use in downstream operations.

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
  aperture = 3L,
  mixed_aperture_level = NULL,
  resround = "nearest"
)
```

## Arguments

- data:

  A data.frame or sf object containing coordinates

- grid:

  A HexGrid object from
  [`hex_grid()`](https://gcol33.github.io/hexify/reference/hex_grid.md).
  If provided, overrides area_km2, resolution, and aperture parameters.

- lon:

  Column name for longitude (ignored if data is sf)

- lat:

  Column name for latitude (ignored if data is sf)

- area_km2:

  Target cell area in km² (mutually exclusive with diagonal). Alias:
  `area` for backwards compatibility.

- diagonal:

  Target cell diagonal (long diagonal) in km

- resolution:

  Grid resolution (0-30). Alternative to area_km2.

- aperture:

  Grid aperture: 3, 4, 7, or "4/3" for mixed (default 3)

- mixed_aperture_level:

  For mixed aperture "4/3": number of aperture-4 levels before switching
  to aperture-3 (default NULL, auto-calculated as resolution/2)

- resround:

  How to round resolution: "nearest", "up", or "down"

## Value

A HexData object containing:

- The input data with cell assignment columns

- The grid specification for downstream operations

Use `as.data.frame(result)` to extract a plain data.frame. Use
`as_sf(result)` to convert to sf object.

Added columns:

- cell_id:

  Stable DGGS cell identifier

- cell_cen_lon:

  Longitude of cell center in degrees

- cell_cen_lat:

  Latitude of cell center in degrees

- cell_area_km2:

  Actual cell area in km²

- cell_diag_km:

  Actual cell long diagonal in km

## Details

For sf objects, coordinates are automatically extracted and transformed
to WGS84 (EPSG:4326) if needed. The geometry column is preserved.

Either `area_km2` (or `area`), `diagonal`, or `resolution` must be
provided unless a `grid` object is supplied.

The HexData return type (default) stores the grid specification so
downstream functions like
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`hexify_cell_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_sf.md),
etc. don't need grid parameters repeated.

## Grid Specification

You can create a grid specification once and reuse it:


    grid <- hex_grid(area_km2 = 1000)
    result1 <- hexify(df1, grid = grid)
    result2 <- hexify(df2, grid = grid)

## See also

[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
grid specification,
[`HexData-class`](https://gcol33.github.io/hexify/reference/HexData-class.md)
for return object details,
[`as_sf`](https://gcol33.github.io/hexify/reference/as_sf.md) for
converting to sf

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
} # }
```
