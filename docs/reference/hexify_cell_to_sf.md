# Convert cell IDs to sf polygons

Creates polygon geometries for hexagonal grid cells from their cell IDs.
Returns an sf object by default, or a data frame for lightweight
workflows.

## Usage

``` r
hexify_cell_to_sf(hex_id, resolution, aperture = 3L, return_sf = TRUE)

hexify_polygons(hex_id, resolution, aperture = 3L, return_sf = TRUE)
```

## Arguments

- hex_id:

  Integer vector of cell identifiers

- resolution:

  Grid resolution level

- aperture:

  Grid aperture: 3, 4, or 7

- return_sf:

  Logical. If TRUE (default), returns sf object with polygon geometries.
  If FALSE, returns data frame with vertex coordinates.

## Value

If return_sf = TRUE: sf object with columns:

- hex_id:

  Cell identifier

- geometry:

  POLYGON geometry (sfc_POLYGON)

If return_sf = FALSE: data frame with columns:

- hex_id:

  Cell identifier

- lon:

  Vertex longitude

- lat:

  Vertex latitude

- order:

  Vertex order (1-7, 7 closes the polygon)

## Details

This function uses a native C++ implementation that is significantly
faster than dggridR's polygon generation, especially for large numbers
of cells.

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)

# Generate some data with hex cells
df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 45))
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# Get polygons as sf object
polys <- hexify_cell_to_sf(result$hex_id, resolution = 10, aperture = 3)

# Plot with sf
library(sf)
plot(st_geometry(polys), col = "lightblue", border = "blue")
} # }
```
