# Convert gridded data or a grid to sf

Methods for
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html).
Gridded data comes back as one feature per row, carrying either the
centre of the cell each row fell in or that cell's boundary. A grid
specification comes back as its global cell set.

## Usage

``` r
# S3 method for class 'HexData'
st_as_sf(
  x,
  ...,
  geometry = c("point", "polygon"),
  shape = c("hexagon", "gosper", "descendants"),
  depth = 3L
)

# S3 method for class 'HexGridInfo'
st_as_sf(x, ...)
```

## Arguments

- x:

  A HexData or HexGridInfo object

- ...:

  For a HexGridInfo, passed on to
  [`grid_global`](https://gillescolling.com/hexify/reference/grid_global.md)

- geometry:

  Type of geometry: "point" (default) or "polygon"

- shape, depth:

  Cell shape for polygon geometry, as for
  [`cell_to_sf`](https://gillescolling.com/hexify/reference/cell_to_sf.md)

## Value

An sf object

## Details

For point geometry, cell centers (cell_cen_lon, cell_cen_lat) are used.
For polygon geometry, cell boundaries are computed using the grid
specification.

## See also

[`cell_to_sf`](https://gillescolling.com/hexify/reference/cell_to_sf.md)
for the cells of named cell IDs,
[`grid_global`](https://gillescolling.com/hexify/reference/grid_global.md)
for every cell of a grid

## Examples

``` r
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

# Get sf points
sf_pts <- st_as_sf(result)

# Get sf polygons
sf_poly <- st_as_sf(result, geometry = "polygon")

# Every cell of a coarse grid
cells <- st_as_sf(hex_grid(resolution = 2))
```
