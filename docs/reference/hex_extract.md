# Extract Raster Values at Hex Cell Centers

Samples raster values at hexagonal cell centers. Faster than
[`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md)
because it only queries cell center points, not full polygons.

## Usage

``` r
hex_extract(raster, grid, cells = NULL, boundary = NULL)
```

## Arguments

- raster:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- grid:

  A HexGridInfo or HexData object specifying the grid.

- cells:

  Optional cell IDs to extract. If `NULL` (default), extracts at all
  cell centers from `grid` (only works for HexData objects or if a
  `boundary` is provided).

- boundary:

  Optional sf polygon to limit extraction extent.

## Value

A data.frame with columns `cell_id`, plus one column per raster layer.

## Details

Requires the `terra` package (in Suggests). The function:

1.  Generates cell centers for the raster extent

2.  Calls `terra::extract(raster, cell_center_matrix)`

3.  Attaches cell IDs

For full zonal statistics (aggregating all pixels within each hex
polygon), use
[`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md)
instead.

## See also

[`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md)
for polygon-based zonal statistics,
[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md) for
creating HexData objects

## Examples

``` r
# \donttest{
if (requireNamespace("terra", quietly = TRUE)) {
  # Create a small synthetic raster
  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = -10, xmax = 10, ymin = 40, ymax = 55)
  terra::values(r) <- runif(100)
  names(r) <- "temperature"

  # Extract at hex cell centers
  g <- hex_grid(area_km2 = 500)
  df <- data.frame(lon = c(0, 5), lat = c(45, 50))
  hd <- hexify(df, lon = "lon", lat = "lat", grid = g)
  hex_extract(r, hd)
}
# }
```
