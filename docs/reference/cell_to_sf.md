# Convert cell IDs to sf polygons

Creates sf polygon geometries for hexagonal grid cells.

## Usage

``` r
cell_to_sf(
  cell_id = NULL,
  grid,
  wrap_dateline = TRUE,
  densify = NULL,
  max_km = NULL,
  shape = c("hexagon", "gosper", "descendants"),
  depth = 3L
)
```

## Arguments

- cell_id:

  Cell IDs (integer64 for ISEA, character for H3). If NULL and x is
  HexData, uses cells from x.

- grid:

  A HexGridInfo or HexData object. If HexData and cell_id is NULL,
  polygons are generated for all cells in the data.

- wrap_dateline:

  Logical. If TRUE (default), calls
  [`sf::st_wrap_dateline()`](https://r-spatial.github.io/sf/reference/st_transform.html)
  to split antimeridian-crossing polygons, which flat maps and planar
  (GEOS) operations need. Set to FALSE for orthographic/globe
  projections where wrapping creates gaps, and for spherical (s2)
  operations: each cell then stays one polygon whose corners lie in
  -180..180, so a cell crossing the antimeridian spans a flat lon/lat
  map.

- densify:

  Points added along each cell edge: every edge is halved until the
  straight line in longitude and latitude between consecutive points
  stays within `densify` times its length of the true edge. `0` keeps
  the corners alone. `NULL` uses 0.001 for ISEA cells and 0 for H3
  cells.

- max_km:

  Longest edge piece, in km on the grid's sphere: edges are halved
  further until no straight piece between consecutive points spans more
  than `max_km`, on top of `densify`. `NULL` (default) sets no limit.
  The bound `densify` sets is relative to each piece's length, so a
  long, gently curved edge can stand as one piece that strays hundreds
  of metres from the true edge; `max_km` bounds that in distance.

- shape:

  `"hexagon"` (default) draws each cell's own boundary. `"gosper"` draws
  it as a Gosper island and `"descendants"` as the outline of its
  descendants; both need an ISEA grid, and `"descendants"` one of
  aperture 7. See Details.

- depth:

  For `shape = "gosper"`, the number of times each edge is replaced,
  multiplying its segments by 3; for `"descendants"`, the number of
  resolutions down, multiplying the cells behind one outline by 7.

## Value

sf object with cell_id and geometry columns

## Details

When called with a HexData object and no cell_id argument, this function
generates polygons for all unique cells in the data, which is useful for
plotting.

An ISEA cell edge is straight on its icosahedron face and curved in
longitude and latitude, so ISEA polygons are densified by default. An H3
cell edge is a great-circle arc between corners, which sf reads exactly
from the corners alone when it uses s2; densify H3 cells for planar
work, or to draw them on a flat map.

A Gosper island replaces every cell edge by three segments
\\1/\sqrt{7}\\ as long, turned by \\\arctan(\sqrt{3}/5)\\, and repeats
this `depth` times; its outline approaches the flowsnake. The
replacement is drawn on the face plane of the projection and moves as
much area out of a cell as into it, so on an equal-area grid an island
has exactly its cell's area, and neighbouring islands share their
boundary, so the islands of a grid tile the sphere. The twelve
pentagonal cells give five-sided islands.

The descendant outline of an aperture-7 cell is the boundary of its
descendants `depth` resolutions down. It covers exactly the cell's area,
the outlines of one resolution tile the sphere, and an outline is the
union of its children's outlines, so they nest across resolutions. The
child lattice of an aperture-7 ISEA grid turns one way at one resolution
and back at the next, so the outline stays close to the hexagon rather
than approaching a Gosper island.

Either shape is drawn around the cell's hexagon: a point near the
boundary can lie in one cell's shape and be assigned to the neighbouring
cell by
[`lonlat_to_cell`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md).

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
grid specifications,
[`st_as_sf`](https://gillescolling.com/hexify/reference/st_as_sf.HexData.md)
for converting HexData to sf

## Examples

``` r
# From grid specification
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
polys <- cell_to_sf(cells, grid)

# From HexData (all cells)
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
polys <- cell_to_sf(grid = result)

# Cells drawn as Gosper islands, and as the outline of their descendants
g7 <- hex_grid(resolution = 3, aperture = 7)
cells <- lonlat_to_cell(c(0, 10), c(45, 50), g7)
islands <- cell_to_sf(cells, g7, shape = "gosper", depth = 3)
outlines <- cell_to_sf(cells, g7, shape = "descendants", depth = 2)
```
