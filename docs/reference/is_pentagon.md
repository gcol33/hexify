# Detect Pentagon Cells

Identifies which cells are pentagons. A hexagonal grid on the
icosahedron has exactly 12 pentagons, one at each icosahedron vertex,
with 5 neighbors instead of 6. A grid on the octahedron has none: its
six vertex cells are squares.

## Usage

``` r
is_pentagon(cell_id, grid)
```

## Arguments

- cell_id:

  Cell IDs to check. integer64 for ISEA, character for H3.

- grid:

  A HexGridInfo or HexData object specifying the grid.

## Value

A logical vector. `TRUE` for pentagon cells, `FALSE` for hexagons.

## Details

**H3 backend:** Uses the vendored H3 `isPentagon` function.

**ISEA backend:** A vertex cell sits at a vertex of the solid, which is
the (i, j) = (0, 0) cell of its quad, and has one side per face meeting
there. Each input cell ID is decoded to its own (quad, i, j) and tested,
rather than by re-deriving each vertex's cell ID (the forward direction
has aperture-specific quirks – e.g. aperture 7's substrate/surrogate
coordinate distinction – that make a single fixed formula for "the (0,0)
cell ID of quad Q" unreliable across apertures).

## See also

[`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md)
for neighbor finding (pentagons have 5 neighbors)

## Examples

``` r
# \donttest{
# H3 pentagon detection
g <- hex_grid(resolution = 1, type = "h3")
cells <- grid_global(g)
pent <- is_pentagon(cells$cell_id, g)
sum(pent)  # Should be 12
# }
```
