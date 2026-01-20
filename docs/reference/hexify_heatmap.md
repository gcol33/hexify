# Create a ggplot2 visualization of hexagonal grid cells

Creates a ggplot2-based visualization of hexagonal grid cells,
optionally colored by a value column. Supports continuous and discrete
color scales, projection transformation, and customizable styling.

## Usage

``` r
hexify_heatmap(
  data,
  value = NULL,
  basemap = NULL,
  crs = NULL,
  colors = NULL,
  breaks = NULL,
  labels = NULL,
  hex_border = "#5D4E37",
  hex_lwd = 0.3,
  hex_alpha = 0.7,
  basemap_fill = "gray90",
  basemap_border = "gray50",
  basemap_lwd = 0.5,
  mask_outside = FALSE,
  aperture = 3L,
  xlim = NULL,
  ylim = NULL,
  title = NULL,
  legend_title = NULL,
  na_color = "gray90",
  theme_void = TRUE
)
```

## Arguments

- data:

  A HexData object from
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md), a
  data frame with cell_id and cell_area columns, or an sf object with
  hexagon polygons.

- value:

  Column name (as string) to use for fill color. If NULL, cells are
  drawn with a uniform fill color. If not specified but data has a
  'count' or 'n' column, that will be used automatically.

- basemap:

  Optional basemap. Can be:

  - `NULL`: No basemap (default)

  - `"world"`: Use built-in `hexify_world` map (low resolution)

  - `"world_hires"`: Use high-resolution map from rnaturalearth
    (requires package)

  - An sf object: User-supplied vector map

- crs:

  Target CRS for the map projection. Can be:

  - A numeric EPSG code (e.g., 4326 for 'WGS84', 3035 for LAEA Europe)

  - A proj4 string

  - An sf crs object

  - NULL to use 'WGS84' (EPSG:4326)

- colors:

  Color palette for the heatmap. Can be:

  - A character vector of colors (for manual scale)

  - A single RColorBrewer palette name (e.g., "YlOrRd", "Greens")

  - NULL to use viridis

- breaks:

  Numeric vector of break points for binning continuous values, or NULL
  for continuous scale. Use `Inf` and `-Inf` for open-ended bins.

- labels:

  Labels for the breaks (length should be one less than breaks). If
  NULL, labels are auto-generated.

- hex_border:

  Border color for hexagons

- hex_lwd:

  Line width for hexagon borders

- hex_alpha:

  Transparency for hexagon fill (0-1)

- basemap_fill:

  Fill color for basemap polygons

- basemap_border:

  Border color for basemap polygons

- basemap_lwd:

  Line width for basemap borders

- mask_outside:

  Logical. If TRUE and basemap is provided, mask hexagon portions that
  fall outside the basemap polygons.

- aperture:

  Grid aperture (default 3), used if data is from hexify()

- xlim:

  Optional x-axis limits (in target CRS units) as c(min, max)

- ylim:

  Optional y-axis limits (in target CRS units) as c(min, max)

- title:

  Plot title

- legend_title:

  Title for the color legend

- na_color:

  Color for NA values

- theme_void:

  Logical. If TRUE (default), use a minimal theme without axes,
  gridlines, or background.

## Value

A ggplot2 object that can be further customized or saved.

## Details

This function provides publication-quality heatmap visualizations of
hexagonal grids using ggplot2. It returns a ggplot object that can be
further customized with standard ggplot2 functions.

## Color Scales

The function supports three types of color scales:

- Continuous:

  Set `breaks = NULL` for a continuous gradient

- Binned:

  Provide `breaks` vector to bin values into categories

- Discrete:

  If `value` column is a factor, discrete colors are used

## Projections

Common projections:

- 4326:

  'WGS84' (unprojected lat/lon)

- 3035:

  LAEA Europe

- 3857:

  Web Mercator

- "+proj=robin":

  Robinson (world maps)

- "+proj=moll":

  Mollweide (equal-area world maps)

## See also

[`plot_grid`](https://gcol33.github.io/hexify/reference/plot_grid.md)
for base R plotting,
[`cell_to_sf`](https://gcol33.github.io/hexify/reference/cell_to_sf.md)
to generate polygons manually

Other visualization:
[`plot_world()`](https://gcol33.github.io/hexify/reference/plot_world.md)

## Examples

``` r
library(hexify)

# Sample data with counts
cities <- data.frame(
  lon = c(16.37, 2.35, -3.70, 12.5, 4.9),
  lat = c(48.21, 48.86, 40.42, 41.9, 52.4),
  count = c(100, 250, 75, 180, 300)
)
result <- hexify(cities, lon = "lon", lat = "lat", area_km2 = 5000)

# Simple plot (uniform fill, no value mapping)
hexify_heatmap(result)

# \donttest{
library(ggplot2)

# With world basemap
hexify_heatmap(result, basemap = "world")

# Heatmap with value mapping
hexify_heatmap(result, value = "count")

# With world basemap and custom colors
hexify_heatmap(result, value = "count",
               basemap = "world",
               colors = "YlOrRd",
               title = "City Density")

# Binned values with custom breaks
hexify_heatmap(result, value = "count",
               basemap = "world",
               breaks = c(-Inf, 100, 200, Inf),
               labels = c("Low", "Medium", "High"),
               colors = c("#fee8c8", "#fdbb84", "#e34a33"))

# Different projection (LAEA Europe)
hexify_heatmap(result, value = "count",
               basemap = "world",
               crs = 3035,
               xlim = c(2500000, 6500000),
               ylim = c(1500000, 5500000))

# Customize further with ggplot2
hexify_heatmap(result, value = "count", basemap = "world") +
  labs(caption = "Data source: Example") +
  theme(legend.position = "bottom")
# }
```
