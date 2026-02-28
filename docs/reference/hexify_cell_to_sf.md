# Convert cell IDs to sf polygons

Creates polygon geometries for hexagonal grid cells from their cell IDs.
Returns an sf object by default, or a data frame for lightweight
workflows.

## Usage

``` r
hexify_cell_to_sf(
  cell_id,
  resolution = NULL,
  aperture = NULL,
  return_sf = TRUE,
  grid = NULL,
  wrap_dateline = TRUE
)
```

## Arguments

- cell_id:

  Integer vector of cell identifiers

- resolution:

  Grid resolution level. Can be omitted if grid is provided.

- aperture:

  Grid aperture: 3, 4, or 7. Can be omitted if grid is provided.

- return_sf:

  Logical. If TRUE (default), returns sf object with polygon geometries.
  If FALSE, returns data frame with vertex coordinates.

- grid:

  Optional HexGridInfo object. If provided, resolution and aperture are
  extracted from it.

- wrap_dateline:

  Logical. If TRUE (default), calls
  [`sf::st_wrap_dateline()`](https://r-spatial.github.io/sf/reference/st_transform.html)
  to split antimeridian-crossing polygons. Set to FALSE for
  orthographic/globe projections where wrapping creates gaps.

## Value

If return_sf = TRUE: sf object with columns:

- cell_id:

  Cell identifier

- geometry:

  POLYGON geometry (sfc_POLYGON)

If return_sf = FALSE: data frame with columns:

- cell_id:

  Cell identifier

- lon:

  Vertex longitude

- lat:

  Vertex latitude

- order:

  Vertex order (1-7, 7 closes the polygon)

## Details

This function uses a native C++ implementation that is significantly
faster than 'dggridR' polygon generation, especially for large numbers
of cells.

For the recommended S4 interface, use
[`cell_to_sf`](https://gcol33.github.io/hexify/reference/cell_to_sf.md)
instead.

## See also

[`cell_to_sf`](https://gcol33.github.io/hexify/reference/cell_to_sf.md)
for the recommended S4 interface

Other sf conversion:
[`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md),
[`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md),
[`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md)

## Examples

``` r
library(hexify)

# Generate some data with hex cells
df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 45))
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# Get polygons as sf object (using HexData)
polys <- cell_to_sf(grid = result)

# Or with explicit parameters
polys <- hexify_cell_to_sf(result@cell_id, resolution = 10, aperture = 3)

# Plot with sf
library(sf)
plot(st_geometry(polys), col = "lightblue", border = "blue")
```
