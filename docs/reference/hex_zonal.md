# Zonal Statistics for Hex Cells

Computes zonal statistics by aggregating all raster pixels falling
within each hexagonal cell polygon. More accurate than
[`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md)
but slower because it requires polygon geometries.

## Usage

``` r
hex_zonal(raster, grid, fun = "mean", boundary = NULL, cells = NULL)
```

## Arguments

- raster:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- grid:

  A HexGridInfo or HexData object specifying the grid.

- fun:

  Summary function name: `"mean"` (default), `"sum"`, `"min"`, `"max"`,
  `"sd"`, or `"count"`.

- boundary:

  Optional sf polygon to limit the analysis extent.

- cells:

  Optional cell IDs. If provided, only these cells are included.

## Value

A data.frame with columns `cell_id`, plus one column per raster layer
containing the aggregated values.

## Details

Requires the `terra` package (in Suggests). The function:

1.  Generates hex polygon geometries via
    [`cell_to_sf()`](https://gillescolling.com/hexify/reference/cell_to_sf.md)

2.  Calls `terra::extract(raster, polygons, fun = fun)`

3.  Joins results back to cell IDs

For point-based extraction (faster, at cell centers only), use
[`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md)
instead.

## See also

[`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md)
for cell-center extraction,
[`cell_to_sf()`](https://gillescolling.com/hexify/reference/cell_to_sf.md)
for hex polygon generation

## Examples

``` r
# \donttest{
if (requireNamespace("terra", quietly = TRUE)) {
  r <- terra::rast(nrows = 100, ncols = 100,
                   xmin = -10, xmax = 10, ymin = 40, ymax = 55)
  terra::values(r) <- runif(10000)
  names(r) <- "temperature"

  g <- hex_grid(area_km2 = 500)
  df <- data.frame(lon = c(0, 5), lat = c(45, 50))
  hd <- hexify(df, lon = "lon", lat = "lat", grid = g)
  hex_zonal(r, hd, fun = "mean")
}
# }
```
