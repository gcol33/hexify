# Plot hexagonal grid cells with optional basemap

Creates a map visualization of hexagonal grid cells. Supports the
built-in world map or user-supplied basemaps (sf vectors or raster
images).

## Usage

``` r
hexify_map(
  data,
  basemap = NULL,
  fill = "steelblue",
  border = "gray30",
  lwd = 1,
  alpha = 0.7,
  basemap_fill = "gray95",
  basemap_border = "gray70",
  basemap_lwd = 1,
  aperture = 3L,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  ...
)
```

## Arguments

- data:

  Data frame from hexify() containing hex_id and hex_area columns, or an
  sf object with hexagon polygons

- basemap:

  Optional basemap. Can be:

  - `NULL`: No basemap (default)

  - `"world"`: Use built-in `hexify_world` map

  - An sf object: User-supplied vector map

  - A SpatRaster (terra) or RasterLayer (raster): User-supplied raster

- fill:

  Fill color for hexagons (single color or column name for mapping)

- border:

  Border color for hexagons

- lwd:

  Line width for hexagon borders

- alpha:

  Transparency for hexagon fill (0-1)

- basemap_fill:

  Fill color for basemap polygons (if vector)

- basemap_border:

  Border color for basemap polygons (if vector)

- basemap_lwd:

  Line width for basemap borders

- aperture:

  Grid aperture (default 3), used if data is from hexify()

- xlim:

  Optional x-axis (longitude) limits as c(min, max)

- ylim:

  Optional y-axis (latitude) limits as c(min, max)

- main:

  Plot title

- ...:

  Additional arguments passed to plot()

## Value

NULL invisibly. Creates a plot as side effect.

## Details

This function provides a simple way to visualize hexagonal grids with
geographic context. For more sophisticated visualizations, use
[`hexify_to_polygons`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md)
to get an sf object and plot with ggplot2, tmap, or other mapping
packages.

The function automatically:

- Converts hexify() output to polygons if needed

- Adjusts aspect ratio for latitude

- Clips basemap to data extent (with buffer)

## Basemap options

- Built-in world map:

  Use `basemap = "world"` for the included simplified Natural Earth map.
  No additional packages required.

- Custom sf vector:

  Pass any sf object as basemap for custom boundaries, regions, or
  detailed coastlines. install.packages("rnaturalearth")

- Raster basemap:

  Pass a SpatRaster (terra) or RasterLayer (raster) for satellite
  imagery or other raster backgrounds. Requires terra or raster package.

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)

# Sample data
cities <- data.frame(
  name = c("Vienna", "Paris", "Madrid"),
  lon = c(16.37, 2.35, -3.70),
  lat = c(48.21, 48.86, 40.42)
)
result <- hexify(cities, lon = "lon", lat = "lat", area = 5000)

# Simple plot without basemap
hexify_map(result)

# With built-in world map
hexify_map(result, basemap = "world")

# Custom colors
hexify_map(result, basemap = "world",
           fill = "steelblue", border = "darkblue",
           basemap_fill = "ivory", basemap_border = "gray50")

# With user-supplied sf basemap
library(rnaturalearth)
europe <- ne_countries(continent = "Europe", returnclass = "sf")
hexify_map(result, basemap = europe)

# Zoom to specific region
hexify_map(result, basemap = "world",
           xlim = c(-10, 25), ylim = c(35, 55))
} # }
```
