# Detect Pentagon Cells

Identifies which cells are pentagons. Any hexagonal tiling of the sphere
must contain exactly 12 pentagons (at the icosahedron vertices).
Pentagon cells have 5 neighbors instead of 6.

## Usage

``` r
is_pentagon(cell_id, grid)
```

## Arguments

- cell_id:

  Cell IDs to check. Numeric for ISEA, character for H3.

- grid:

  A HexGridInfo or HexData object specifying the grid.

## Value

A logical vector. `TRUE` for pentagon cells, `FALSE` for hexagons.

## Details

**H3 backend:** Uses the vendored H3 `isPentagon` function.

**ISEA backend:** The 12 pentagons are located at icosahedron vertices,
which are always the (i, j) = (0, 0) cell of their quad. Pentagon status
is checked by decoding each input cell ID's own (quad, i, j) coordinates
via
[`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)'s
underlying grid math and testing whether i and j are both zero, rather
than by re-deriving each vertex's cell ID (the forward direction has
aperture-specific quirks – e.g. aperture 7's substrate/surrogate
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
