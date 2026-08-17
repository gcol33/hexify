# Grid Distance Between Cells

Computes the grid distance (minimum number of hops) between pairs of
hexagonal cells. This is a discrete distance measured in cell steps, not
a geodesic distance.

## Usage

``` r
hex_distance(cell_a, cell_b, grid)
```

## Arguments

- cell_a, cell_b:

  Cell IDs. Must be the same length (pairs) or one of them length 1
  (broadcast). For H3 grids, character vectors. For ISEA grids, numeric
  vectors.

- grid:

  A HexGridInfo or HexData object specifying the grid.

## Value

An integer vector of grid distances. `NA` where the distance cannot be
computed (e.g., pentagon path issues in H3).

## Details

**H3 backend:** Uses the vendored H3 `gridDistance` function. Returns
the exact shortest path length in cell hops.

**ISEA backend:** For cells on the same quad, computes the
cube-coordinate distance: `max(|di|, |dj|, |di + dj|)`. Cross-quad
distances use BFS expansion and may return `NA` for very distant cells.

For geodesic (geographic) distances between cell centers, convert to
lon/lat with
[`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
and use
[`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html).

## See also

[`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md)
for finding cells within a given distance,
[`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
for geographic coordinates

## Examples

``` r
# \donttest{
# H3 grid distance
g <- hex_grid(resolution = 5, type = "h3")
a <- lonlat_to_cell(10, 50, g)
b <- lonlat_to_cell(10.1, 50.1, g)
hex_distance(a, b, g)
# }
```
