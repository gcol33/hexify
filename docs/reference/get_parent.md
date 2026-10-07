# Get parent cell

Returns the parent cell at a coarser resolution.

## Usage

``` r
get_parent(cell_id, grid, levels = 1L, overlapping = FALSE)
```

## Arguments

- cell_id:

  Cell IDs: integer64 for ISEA grids (whole numbers below 2^53 and digit
  strings are accepted), character for H3 grids

- grid:

  A HexGridInfo or HexData object

- levels:

  Number of levels up (default 1)

- overlapping:

  If `TRUE`, return every cell at the coarser resolution that overlaps
  the cell, not only the one holding its centre.

## Value

Parent cell IDs, of the type the grid's cell IDs have. With
`overlapping = TRUE`, a list with one vector per input cell, the parent
holding the cell's centre first.

## Details

A cell's parent is the coarser cell holding its centre. The hexagons of
successive resolutions do not nest, so a cell can also reach into one or
two neighbours of that parent: an aperture-3 cell centred on a parent
corner lies a third in each of three parents, an aperture-4 cell centred
on a parent edge half in each of two, and an aperture-7 cell off the
parent's centre 11/12 in its parent and 1/12 in one neighbour.
`overlapping = TRUE` returns all of them, as the `parent()` query of OGC
Topic 21 does when `inheritID` is false; the default is its
`inheritID = true` answer.

Every part of a cell that lies in a coarser cell holds a corner of the
cell, so the overlapping cells are the coarser cells holding the cell's
corners, each taken a small step towards the centre.

## Examples

``` r
grid <- hex_grid(resolution = 10)
child_cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
parent_cells <- get_parent(child_cells, grid)

# Every coarser cell a cell reaches into
g3 <- hex_grid(resolution = 5, aperture = 3)
cells <- lonlat_to_cell(c(0, 10), c(45, 50), g3)
get_parent(cells, g3, overlapping = TRUE)
```
