# Ray glyphs for cell values

Draws one ray per cell from its centre whose angle shows a value, as in
Carr et al. (1997): the ray points down for the smallest value, level
for the middle of the range and up for the largest. A confidence arc at
the ray's tip spans the angles of a lower and an upper bound. A second
variable gets a ray of its own pointing to the left from the same
centre, so one glyph shows two values.

## Usage

``` r
hex_rays(
  cell_id,
  value,
  grid,
  lower = NULL,
  upper = NULL,
  value2 = NULL,
  lower2 = NULL,
  upper2 = NULL,
  limits = NULL,
  limits2 = NULL,
  length = 0.85
)
```

## Arguments

- cell_id:

  Cell IDs, each listed once.

- value:

  Values shown by the right-hand rays.

- grid:

  A HexGridInfo or HexData object.

- lower, upper:

  Optional bounds of `value` (such as a confidence interval), drawn as
  an arc at the tip of the ray.

- value2, lower2, upper2:

  Optional second variable and its bounds, drawn to the left.

- limits:

  Values drawn straight down and straight up: `c(min, max)`. `NULL` uses
  the range of the values and bounds; with a second variable both share
  it unless `limits2` is given.

- limits2:

  Limits of the second variable.

- length:

  Length of a ray as a fraction of the distance from the cell centre to
  the middle of a wall.

## Value

An sf object of LINESTRINGs with columns `cell_id`, `side` (`"right"` or
`"left"`) and `part` (`"ray"` or `"arc"`).

## Details

Rays run in the plane tangent to the sphere at the cell centre, with up
towards north, and are returned as lines in longitude and latitude,
ready for any map of the cells.

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`hex_triangles`](https://gillescolling.com/hexify/reference/hex_triangles.md)
for the change towards each neighbour

## Examples

``` r
g <- hex_grid(resolution = 4, aperture = 3)
cells <- lonlat_to_cell(c(5, 10, 15), c(45, 47, 49), g)
trend <- c(-1, 0.2, 1.4)
rays <- hex_rays(cells, trend, g, lower = trend - 0.5, upper = trend + 0.5)
plot(sf::st_geometry(cell_to_sf(cells, g)))
plot(sf::st_geometry(rays), add = TRUE, col = "#D55E00")
```
