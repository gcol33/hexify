# Create a ggplot2 visualization of HexData

Generates a ggplot2 object for HexData, supporting fill mapping,
basemaps, and advanced customization.

## Usage

``` r
autoplot.HexData(
  object,
  basemap = TRUE,
  basemap_fill = "gray90",
  basemap_border = "gray50",
  basemap_lwd = 0.3,
  grid_border = "#5D4E37",
  grid_lwd = 0.4,
  grid_alpha = 0.7,
  fill = NULL,
  show_points = FALSE,
  point_size = 1,
  point_color = "red",
  crop = TRUE,
  crop_expand = 0.1,
  title = NULL,
  legend_title = NULL,
  ...
)
```

## Arguments

- object:

  A HexData object from
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)

- basemap:

  Basemap specification (see `plot.HexData`)

- basemap_fill:

  Fill color for basemap

- basemap_border:

  Border color for basemap

- basemap_lwd:

  Line width for basemap

- grid_border:

  Border color for grid cells

- grid_lwd:

  Line width for cell borders

- grid_alpha:

  Transparency for cell fill

- fill:

  Column name for fill mapping (optional)

- show_points:

  Show original points

- point_size:

  Size of points

- point_color:

  Color of points

- crop:

  Crop to data extent

- crop_expand:

  Expansion factor for crop

- title:

  Plot title

- legend_title:

  Legend title

- ...:

  Additional arguments (ignored)

## Value

A ggplot object that can be further customized

## Details

Requires ggplot2 package. The returned object can be modified using
standard ggplot2 functions like `+`,
[`theme()`](https://ggplot2.tidyverse.org/reference/theme.html), etc.

## See also

[`plot,HexData,missing-method`](https://gcol33.github.io/hexify/reference/plot-HexData-missing-method.md)
for base R plotting

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Basic plot
autoplot(result)

# With fill mapping
result$data$count <- sample(1:100, nrow(result$data))
autoplot(result, fill = "count") +
  scale_fill_viridis_c()

# Customize
autoplot(result, basemap = TRUE) +
  theme_minimal() +
  labs(title = "Hexified Data")
} # }
```
