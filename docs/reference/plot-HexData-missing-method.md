# Plot HexData objects

Default plot method for HexData objects. Draws hexagonal grid cells with
an optional basemap.

## Usage

``` r
# S4 method for class 'HexData,missing'
plot(
  x,
  y,
  basemap = TRUE,
  basemap_fill = "ivory",
  basemap_border = "gray50",
  basemap_lwd = 0.5,
  grid_fill = "#3B9AB2",
  grid_border = "#1A1A1A",
  grid_lwd = 0.8,
  grid_alpha = 0.6,
  fill = NULL,
  show_points = FALSE,
  point_size = 1,
  point_color = "red",
  crop = TRUE,
  crop_expand = 0.1,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A HexData object from
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)

- y:

  Ignored (for S4 method compatibility)

- basemap:

  Basemap specification:

  - `TRUE` or `"world"`: Use built-in world map

  - `FALSE` or `NULL`: No basemap

  - sf object: Custom basemap

- basemap_fill:

  Fill color for basemap (default "ivory")

- basemap_border:

  Border color for basemap (default "gray50")

- basemap_lwd:

  Line width for basemap borders (default 0.5)

- grid_fill:

  Fill color for grid cells (default "#3B9AB2" - teal)

- grid_border:

  Border color for grid cells (default "#1A1A1A" - near-black)

- grid_lwd:

  Line width for cell borders (default 0.8)

- grid_alpha:

  Transparency for cell fill (0-1, default 0.6)

- fill:

  Column name for fill mapping (optional)

- show_points:

  Show original points on top of cells (default FALSE)

- point_size:

  Size of points if shown (default 1)

- point_color:

  Color of points (default "red")

- crop:

  Crop to data extent (default TRUE)

- crop_expand:

  Expansion factor for crop (default 0.1)

- main:

  Plot title

- ...:

  Additional arguments passed to base plot()

## Value

Invisibly returns the HexData object

## Details

This function generates polygon geometries for the cells present in the
data and plots them. Polygons are computed on demand, not stored, to
minimize memory usage.

## See also

[`autoplot.HexData`](https://gcol33.github.io/hexify/reference/autoplot.HexData.md)
for ggplot2 plotting

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(lon = runif(100, -10, 10), lat = runif(100, 40, 50))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Basic plot
plot(result)

# With world basemap
plot(result, basemap = TRUE)

# Custom styling
plot(result, basemap = TRUE,
     grid_fill = "lightblue", grid_border = "darkblue",
     basemap_fill = "ivory")

# Show original points
plot(result, basemap = TRUE, show_points = TRUE)
} # }
```
