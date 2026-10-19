# HexData Class

An S4 class representing hexified data. Contains the original user data
plus cell assignments from the hexification process.

## Usage

``` r
# S4 method for class 'HexData'
grid_info(x)

# S4 method for class 'HexData'
cells(x)

# S4 method for class 'HexData'
n_cells(x)

# S4 method for class 'HexData'
nrow(x)

# S4 method for class 'HexData'
ncol(x)

# S4 method for class 'HexData'
dim(x)

# S4 method for class 'HexData'
names(x)

# S4 method for class 'HexData'
x$name

# S4 method for class 'HexData'
x$name <- value

# S4 method for class 'HexData'
x[i, j, ..., drop = FALSE]

# S4 method for class 'HexData'
x[[i]]

# S4 method for class 'HexData,ANY,missing'
x[[i, j]] <- value

# S4 method for class 'HexData'
show(object)

# S4 method for class 'HexData'
summary(object, ...)

# S3 method for class 'hexify_data_summary'
print(x, ...)

# S4 method for class 'HexData'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S4 method for class 'HexData'
as.list(x, ...)

# S3 method for class 'HexData'
as_tibble(x, ...)

# S3 method for class 'HexData'
st_as_sf(
  x,
  ...,
  geometry = c("point", "polygon"),
  shape = c("hexagon", "gosper", "descendants"),
  depth = 3L
)

# S3 method for class 'HexData'
st_crs(x, ...)

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
  point_alpha = 1,
  jitter = TRUE,
  crop = TRUE,
  crop_expand = 0.1,
  main = NULL,
  shape = c("hexagon", "gosper", "descendants"),
  depth = 3L,
  ...
)
```

## Arguments

- x:

  A HexData object; for [`print()`](https://rdrr.io/r/base/print.html),
  a data summary as [`summary()`](https://rdrr.io/r/base/summary.html)
  returns one

- name:

  Column name

- value:

  Replacement value

- i, j:

  Row/column indices

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), passed
  on to base [`plot()`](https://rdrr.io/r/graphics/plot.default.html);
  ignored by the other methods

- drop:

  Logical, whether to drop dimensions

- object:

  A HexData object

- row.names:

  Optional row names

- optional:

  Logical (ignored)

- geometry:

  For
  [`st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html):
  "point" (default) gives each row the centre of its cell, "polygon" the
  cell's boundary

- shape, depth:

  Cell shape for polygon geometry and for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), as for
  [`cell_to_sf`](https://gillescolling.com/hexify/reference/cell_to_sf.md)

- y:

  Ignored (for S4 method compatibility)

- basemap:

  Basemap specification:

  - `TRUE` or `"world"`: Use built-in world map

  - `FALSE` or `NULL`: No basemap

  - `"world_hires"`: High-resolution map from 'rnaturalearth' (requires
    the package)

  - An sf object: User-supplied vector map

  - A SpatRaster: User-supplied raster, drawn in grayscale

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

- point_alpha:

  Transparency for points (0-1, default 1)

- jitter:

  If TRUE (default), points are randomly scattered within their assigned
  hexagon. If FALSE, all points in a cell are plotted at the cell
  centroid instead.

- crop:

  Crop to data extent (default TRUE)

- crop_expand:

  Expansion factor for crop (default 0.1)

- main:

  Plot title

## Value

- [`show()`](https://rdrr.io/r/methods/show.html): the object,
  invisibly, after printing two lines: rows, columns and cells, and the
  grid

- [`summary()`](https://rdrr.io/r/base/summary.html): a list of class
  `hexify_data_summary` carrying `rows`, `columns`, `column_names`,
  `n_cells`, `type`, the `grid` summary and a `preview` of the first
  rows; printing it lists all of these and returns it invisibly

- [`grid_info()`](https://gillescolling.com/hexify/reference/grid_info.md):
  HexGridInfo object containing grid specification

- [`cells()`](https://gillescolling.com/hexify/reference/cells.md):
  vector of unique cell IDs

- [`n_cells()`](https://gillescolling.com/hexify/reference/n_cells.md):
  integer count of unique cells

- [`nrow()`](https://rdrr.io/r/base/nrow.html),
  [`ncol()`](https://rdrr.io/r/base/nrow.html),
  [`dim()`](https://rdrr.io/r/base/dim.html): integer dimensions

- [`names()`](https://rdrr.io/r/base/names.html): character vector of
  column names (including virtual cell columns)

- `$`, `[[`: the requested column or cell data as a vector

- `$<-`, `[[<-`: the modified HexData object

- `[`: subsetted HexData object or extracted data

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): data
  frame with original data plus cell columns

- [`as.list()`](https://rdrr.io/r/base/list.html): named list containing
  data, grid, cell_id, and cell_center

- `as_tibble()`: the original data (without geometry) as a tibble

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the object,
  invisibly

- [`st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html):
  an sf object, one feature per row

- [`st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html):
  the grid's coordinate reference system, as for a HexGridInfo

## Details

HexData objects are created by
[`hexify`](https://gillescolling.com/hexify/reference/hexify.md). The
original data is preserved in the `data` slot, while cell assignments
are stored separately in `cell_id` and `cell_center`.

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
get a combined data frame with cell columns.

Unlike `[.data.frame`, `drop` defaults to `FALSE`: selecting a single
column returns a `HexData` object (preserving `grid`/
`cell_id`/`cell_center`) rather than dropping to a bare vector. Pass
`drop = TRUE` explicitly to get data.frame-style dropping.

## Slots

- `data`:

  Data frame or sf object. The original user data (untouched).

- `grid`:

  HexGridInfo object. The grid specification used.

- `cell_id`:

  Cell IDs for each row of data:
  [`bit64::integer64`](https://bit64.r-lib.org/reference/bit64-package.html)
  for ISEA grids, character for H3 grids.

- `cell_center`:

  Matrix. Two-column matrix (lon, lat) of cell centers.

## Plotting data

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the cells
the data falls in, with an optional basemap. It builds polygons for the
cells present in the data when it draws them, rather than storing them,
to keep memory use low.

## See also

[`hexify`](https://gillescolling.com/hexify/reference/hexify.md) for
creating HexData objects,
[`HexGridInfo-class`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
for grid specifications

[`hexify_heatmap`](https://gillescolling.com/hexify/reference/hexify_heatmap.md)
for ggplot2 plotting

## Examples

``` r
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55), n = 1:3)
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
result
summary(result)
summary(result)$column_names
cells(result)
n_cells(result)
head(as.data.frame(result))
result$cell_id
result[1:2, ]
st_crs(result)

# sf points at the cell centres, or the cells as polygons
sf_pts <- st_as_sf(result)
sf_poly <- st_as_sf(result, geometry = "polygon")

if (requireNamespace("tibble", quietly = TRUE)) {
  tibble::as_tibble(result)
}
# \donttest{
df <- data.frame(lon = runif(50, -5, 5), lat = runif(50, 45, 50))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)

# Basic plot
plot(result, basemap = FALSE)

# With basemap and custom styling
plot(result, grid_fill = "lightblue", grid_border = "darkblue")
# }
```
