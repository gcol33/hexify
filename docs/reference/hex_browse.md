# Interactive Hex Map

Opens an interactive leaflet map with hexagonal cell polygons. Cells can
be colored by a data column, with popups showing cell information on
click.

## Usage

``` r
hex_browse(
  hex_data,
  grid = NULL,
  value = NULL,
  palette = "viridis",
  opacity = 0.7
)
```

## Arguments

- hex_data:

  A HexData object, or a data.frame with a `cell_id` column.

- grid:

  A HexGridInfo object. Required if `hex_data` is a data.frame.

- value:

  Optional column name (character) to color cells by. If `NULL`, cells
  are colored uniformly.

- palette:

  Color palette name for continuous values. Default `"viridis"`.

- opacity:

  Fill opacity (0-1). Default 0.7.

## Value

A `leaflet` map object (can be printed or embedded in Shiny).

## Details

Requires the `leaflet` package (in Suggests). The map is built using
[`as_sf()`](https://gillescolling.com/hexify/reference/as_sf.md) to
generate polygon geometries, then rendered as a leaflet choropleth.

## See also

[`hexify_heatmap()`](https://gillescolling.com/hexify/reference/hexify_heatmap.md)
for static ggplot2 maps,
[`as_sf()`](https://gillescolling.com/hexify/reference/as_sf.md) for sf
conversion

## Examples

``` r
# \donttest{
if (requireNamespace("leaflet", quietly = TRUE)) {
  df <- data.frame(
    lon = runif(50, -5, 5),
    lat = runif(50, 45, 55),
    value = rnorm(50)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 500)
  hex_browse(hd, value = "value")
}
# }
```
