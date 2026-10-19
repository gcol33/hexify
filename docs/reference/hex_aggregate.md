# Aggregate Cell Values to a Coarser Resolution

Adds up, or averages, values given on cells into the cells `levels`
resolutions coarser. The hexagons of successive resolutions do not nest,
so a cell can lie partly in two or three coarser cells; `rule` says how
its value is divided among them.

## Usage

``` r
hex_aggregate(
  cell_id,
  value,
  grid,
  levels = 1L,
  rule = c("area", "centre"),
  measure = c("total", "mean"),
  na.rm = FALSE
)
```

## Arguments

- cell_id:

  Cell IDs, each listed once.

- value:

  Numeric vector with one value per cell, or a numeric matrix or data
  frame with one row per cell, each column aggregated on its own.

- grid:

  A HexGridInfo or HexData object at the resolution of `cell_id`.

- levels:

  Number of resolutions up.

- rule:

  `"area"` divides each cell among the coarser cells it overlaps by the
  share of its area in each; `"centre"` gives each cell wholly to the
  coarser cell holding its centre, the one
  [`get_parent`](https://gillescolling.com/hexify/reference/get_parent.md)
  returns.

- measure:

  `"total"` adds the values up, for counts and other amounts; `"mean"`
  takes their mean weighted by the area each cell contributes, for
  densities and other values per unit area.

- na.rm:

  If `TRUE`, a cell whose value is `NA` counts as absent from that
  column.

## Value

A data frame with one row per coarser cell that receives a share of any
listed cell, sorted by `cell_id`, and the aggregated values: a column
`value` for a vector, the columns of `value` otherwise. With
`measure = "mean"`, a coarser cell whose contributing cells are all `NA`
gets `NA`.

## Details

**The two rules.** Write \\\|c \cap p\| / \|c\|\\ for the share of cell
\\c\\'s area that lies in the coarser cell \\p\\. The area rule gives
\\p\\ that share of \\c\\'s value, which assumes the amount is spread
evenly over \\c\\; the centre rule gives \\p\\ all of it when \\p\\
holds \\c\\'s centre and nothing otherwise. A cell's shares add up to
one under both, so both keep the total of the values. Only the area rule
keeps where the amount lies: a uniform density gives every hexagonal
coarser cell the same total under it, and each vertex cell (the twelve
pentagons of the icosahedron) its sides over six of that. A mean is the
total of value times area over the total area, both aggregated by the
same rule.

**Shares on the lattice.** Within a face, the cells of one resolution
are the Voronoi cells of a hexagonal lattice and those of the next
coarser resolution the Voronoi cells of a sublattice, of index 3, 4 or
7, the aperture of the step. After an aperture-3 step the cells are
centred on the parent centres and corners: the three parent edges
meeting at a corner run through every second edge midpoint of the cell
centred there and cut it into thirds. After an aperture-4 step the cells
are centred on the parent centres and edge midpoints: a parent edge runs
through the centre and two opposite corners of the cell on its midpoint
and cuts it in half. After an aperture-7 step every centre lies inside a
parent. Take the finer lattice as the Eisenstein integers
\\Z\[\omega\]\\, \\\omega = e^{2\pi i/3}\\, with unit spacing, so a cell
has circumradius \\s = 1/\sqrt{3}\\. The parent lattice is \\\mu
Z\[\omega\]\\ with \\\mu = 3 + \omega\\, \\\|\mu\|^2 = 7\\, turned by
\\\arg \mu = \arctan(\sqrt{3}/5)\\ (the mirror image \\3 + \bar\omega\\
at the other twist). The centre child, of circumradius \\s\\, lies
inside its parent, of inradius \\\sqrt{7}/2\\. The parent edge facing
the neighbour at \\\mu\\ runs from the parent corner \\\mu(1 - \omega)/3
= 1 + (1 - \omega)/3\\, which is a corner \\A\\ of the ring child
centred at 1, through the parent's edge midpoint \\\mu/2 = (1 + (2 +
\omega))/2\\, which is the midpoint of the edge the children at 1 and
\\2 + \omega\\ share. That edge is the second from \\A\\, so the part of
the child at 1 beyond its parent is the triangle of \\A\\, the next
corner \\B\\ and \\\mu/2\\: sides \\s\\ and \\s/2\\ at 120 degrees, area
\\\sqrt{3}s^2/8\\ against \\3\sqrt{3}s^2/2\\ for the cell, one twelfth.
The half-turn about \\\mu/2\\ swaps the two parents and the two
children, so the edge cuts one twelfth off the child at \\2 + \omega\\
the other way. Each of the six ring children lies 11/12 in its parent
and 1/12 in one neighbour, and the parent's area checks: \\1 + 6 \cdot
11/12 + 6 \cdot 1/12 = 7\\.

A Hex9 step (aperture 9 on the octahedron, Griffin 2026) divides each
triangle of the solid's lattice into nine. A cell is the six lattice
triangles around its centre, so a parent's edges are lines of the
coarser lattice and therefore of the finer one. Six of a parent's nine
children lie inside it; the other three are centred on its edges, and
the edge through such a child's centre cuts it into its two
half-hexagons, three triangles each: half in each of two parents. The
parent's area checks: \\6 + 6 \cdot 1/2 = 9\\ children's worth.

The vertices of the solid are cell centres at every resolution, and the
turn by the angle deficit about a vertex maps both lattices onto
themselves, so the cells around a vertex are cut as in the plane; the
cell at a vertex lies wholly inside the vertex cell above it. On Hex9
the vertices are lattice points no cell is centred on, and the cells
there are cut along lattice lines as everywhere. Snyder's projection and
the vertex-oriented one (`"ivea"`) are equal-area, so these plane shares
are the shares on the sphere, for every cell, across quad and face edges
and at the vertices.
[`get_parent`](https://gillescolling.com/hexify/reference/get_parent.md)
with `overlapping = TRUE` finds the coarser cells a cell overlaps, its
own parent first, and the share follows from their number and the
aperture of the step: one parent, 1; three after an aperture-3 step, 1/3
each; two after aperture 4, 1/2 each; two after aperture 7, 11/12 to the
parent and 1/12 to the neighbour; two after a Hex9 step, 1/2 each, the
parent holding the child's mode-0 half first.

**Several levels.** `levels > 1` applies the one-level rule level by
level, each step with the shares of its own aperture: a grid built with
`aperture = c(3, 7, 4)` aggregates from resolution 3 with halves, then
twelfths, then thirds. A value at an intermediate resolution is again
taken as spread evenly over its cell, so several levels up the area rule
is the one-level rule repeated rather than the share of a fine cell's
area in a coarse cell.

**Fuller's projection and H3.** Neither is equal-area, so the plane
shares do not carry over to the sphere, and each share is measured. Both
cells are clipped against each other where they are straight, which
makes the pieces exact: an ISEA-family cell on each face plane it
covers, an H3 cell on the gnomonic plane at the finer cell's centre,
where its great-circle edges are straight. An H3 piece is a spherical
polygon whose area is a sum of spherical triangles, exact to rounding. A
face-plane piece has curved edges on the sphere; its area is that of the
polygon through equally spaced points of its edges, extrapolated to
infinitely many points (Romberg), which leaves an error near rounding.
The share is a piece's area over that of all the cell's pieces. Measured
this way on Snyder's projection, the shares come out as the lattice's.

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`get_parent`](https://gillescolling.com/hexify/reference/get_parent.md)
for the coarser cells a cell overlaps,
[`hex_summarize`](https://gillescolling.com/hexify/reference/hex_summarize.md)
to aggregate points into cells

## Examples

``` r
g <- hex_grid(resolution = 4, aperture = 3)
ids <- grid_global(g)$cell_id
set.seed(1)
counts <- rpois(length(ids), 5)

area <- hex_aggregate(ids, counts, g)
centre <- hex_aggregate(ids, counts, g, rule = "centre")
c(sum(counts), sum(area$value), sum(centre$value))

# A uniform density: every parent gets its own area, the hexagons one
# value and the twelve pentagons 5/6 of it
density <- cell_area(ids, g)
range(hex_aggregate(ids, density, g, levels = 2)$value)

# Values per unit area average instead
hex_aggregate(ids, data.frame(a = counts, b = 2 * counts), g,
              measure = "mean")[1:3, ]
```
