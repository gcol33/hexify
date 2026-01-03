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
  clip_basemap = TRUE,
  basemap_fill = "gray90",
  basemap_border = "gray50",
  basemap_lwd = 0.5,
  grid_fill = "#E69F00",
  grid_border = "#5D4E37",
  grid_lwd = 0.8,
  grid_alpha = 0.7,
  fill = NULL,
  show_points = FALSE,
  point_size = "auto",
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

- clip_basemap:

  Clip basemap to data extent (default TRUE). Clipping temporarily
  disables S2 spherical geometry to avoid edge-crossing errors.

- basemap_fill:

  Fill color for basemap (default "gray90")

- basemap_border:

  Border color for basemap (default "gray50")

- basemap_lwd:

  Line width for basemap borders (default 0.5)

- grid_fill:

  Fill color for grid cells (default "#E69F00" - amber/orange)

- grid_border:

  Border color for grid cells (default "#5D4E37" - dark brown)

- grid_lwd:

  Line width for cell borders (default 0.8)

- grid_alpha:

  Transparency for cell fill (0-1, default 0.7)

- fill:

  Column name for fill mapping (optional)

- show_points:

  Show original points on top of cells (default FALSE). Points are
  jittered within their assigned hexagon.

- point_size:

  Size of points. Can be:

  - A number (direct cex value)

  - A preset defining what fraction of a hex cell one point covers:
    "tiny" (~2\\ "large" (~20\\

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

[`hexify_heatmap`](https://gcol33.github.io/hexify/reference/hexify_heatmap.md)
for ggplot2 plotting

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(lon = runif(100, -10, 10), lat = runif(100, 40, 50))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Basic plot (basemap shown by default)
plot(result)

# Without basemap
plot(result, basemap = FALSE)

# Custom styling
plot(result,
     grid_fill = "lightblue", grid_border = "darkblue",
     basemap_fill = "ivory")

# Show jittered points (auto-sized based on density)
plot(result, show_points = TRUE)

# Control point size with presets
plot(result, show_points = TRUE, point_size = "small")
plot(result, show_points = TRUE, point_size = "large")
} # }
```
