# The parts of cells on a net layout

Cuts each cell into its parts on the pieces of a net layout and places
them in the plane. A cell on one piece has one part. A cell that crosses
from one piece to another stays whole on the map where the two pieces
are joined, and is split where a seam of the net runs through it: its
parts then lie apart, and Carr et al. (1997) call it a cloned cell. A
cell on a piece the layout shows twice also has a part on each copy.

## Usage

``` r
net_cells(layout, grid, cells = NULL)
```

## Arguments

- layout:

  A `hexify_net` object from
  [`net_layout`](https://gillescolling.com/hexify/reference/net_layout.md).

- grid:

  The HexGridInfo object the layout was built for.

- cells:

  Cell IDs; `NULL` for every cell of the grid.

## Value

A data frame with one row per part: `cell_id`, `piece` (position in
`layout$pieces`), `group` (parts of a cell that the map joins share a
group, numbered from 1 within the cell), `area` (map area in units of a
face edge squared) and the centroid `x`, `y`. A cell whose parts fall
into more than one group is cut by a seam.

## Details

Cell walls are straight on the faces of the solid and the layout moves
pieces rigidly (or, for `"rhombic"`, by one shear that keeps areas to
one constant), so a part's map area is exact. Snyder's projection is
equal-area, so on a grid built on it the map areas of a cell's parts add
up to the cell's share of the sphere: a face of the net has area
\\\sqrt{3}/4\\, and \\4\pi/n\\ of the unit sphere for a solid of \\n\\
faces.

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`net_layout`](https://gillescolling.com/hexify/reference/net_layout.md),
[`net_project`](https://gillescolling.com/hexify/reference/net_project.md)

## Examples

``` r
g <- hex_grid(resolution = 3, aperture = 3)
net <- net_layout(g, "plane")
parts <- net_cells(net, g)
cut <- tapply(parts$group, parts$cell_id, max) > 1
sum(cut)

# On Snyder's equal-area projection every hexagon has the same map area,
# however the net cuts it
area <- tapply(parts$area, parts$cell_id, sum)
hexagon <- !is_pentagon(as.numeric(names(area)), g)
range(area[hexagon])
```
