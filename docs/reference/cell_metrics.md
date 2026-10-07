# Shape metrics of cells: area, perimeter and compactness

Measures each cell on the sphere: its area, the length of its boundary
and how close its shape is to a circle.

## Usage

``` r
cell_metrics(cell_id = NULL, grid)
```

## Arguments

- cell_id:

  Cell IDs. When `NULL`, the cells of a HexData `grid`, or every cell of
  a HexGridInfo `grid`.

- grid:

  A HexGridInfo or HexData object.

## Value

Data frame with one row per `cell_id`: `cell_id`, `area_km2` (as
[`cell_area`](https://gillescolling.com/hexify/reference/cell_area.md)
gives it), `normalized_area`, `perimeter_km`, `compactness` and `ipq`.

## Details

The perimeter follows each wall as the grid draws it. An ISEA wall is
straight on its face of the solid and curved on the sphere, so it is
followed on the sphere until its length is exact to about 1e-8. An H3
wall is a great-circle arc between corners, with an extra corner where
it crosses an edge of the icosahedron.

Compactness (White et al. 1998) is the perimeter of a spherical cap with
the cell's area over the cell's perimeter, \\\sqrt{4\pi a - a^2/r^2} /
p\\ for area \\a\\, perimeter \\p\\ and radius \\r\\. A circle scores 1;
a small regular hexagon 0.952, pentagon 0.930, square 0.886 and triangle
0.778.

The normalized area is the cell's area over the mean area of the grid's
cells, the body's surface over the cell count; an equal-area grid of
\\N\\ cells gives \\N / (N - 2)\\ for its hexagons, 5/6 of that for the
pentagons of the icosahedron and 4/6 for the squares of the octahedron.
The isoperimetric quotient `ipq` is \\4\pi a / p^2\\, for a circle 1, a
small regular hexagon 0.907, square 0.785 and triangle 0.605. Kmoch et
al. (2022) compare grids by these two measures, reading area and
perimeter in an equal-area plane centred on each cell; here both are
read on the sphere, which adds \\(a/(rp))^2\\ to `compactness`^2, a few
parts in a million for cells of a few thousand square kilometres.

## References

White, D., Kimerling, A. J., Sahr, K. and Song, L. (1998). Comparing
area and shape distortion on polyhedral-based recursive partitions of
the sphere. International Journal of Geographical Information Science
12(8), 805-827.
[doi:10.1080/136588198241518](https://doi.org/10.1080/136588198241518)

Kmoch, A., Vasilyev, I., Virro, H. and Uuemaa, E. (2022). Area and shape
distortions in open-source discrete global grid systems. Big Earth Data
6(3), 256-275.
[doi:10.1080/20964471.2022.2094926](https://doi.org/10.1080/20964471.2022.2094926)

## See also

[`wall_metrics`](https://gillescolling.com/hexify/reference/wall_metrics.md)
for the spacing of neighbouring cells,
[`cell_area`](https://gillescolling.com/hexify/reference/cell_area.md)

## Examples

``` r
g <- hex_grid(resolution = 3, aperture = 3)
m <- cell_metrics(grid = g)
summary(m$compactness)

# Fuller's projection is not equal-area: compare the spread of areas
gf <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
mf <- cell_metrics(grid = gf)
sd(mf$area_km2) / mean(mf$area_km2)
```
