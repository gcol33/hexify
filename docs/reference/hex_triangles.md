# The triangles of cells towards each neighbour

Cuts each cell into one triangle per wall, from the centre to that wall:
six for a hexagon, five for a pentagon. With values, each triangle
carries the change from its cell's value to the neighbour's across the
wall, so colouring the triangles shows the change from cell to cell in
every direction (Carr et al. 1997, after K. and R. Keister).

## Usage

``` r
hex_triangles(cell_id, grid, value = NULL)
```

## Arguments

- cell_id:

  Cell IDs, each listed once.

- grid:

  A HexGridInfo or HexData object.

- value:

  Optional values, one per cell. A neighbour not in `cell_id` or with an
  `NA` value gives an `NA` change.

## Value

An sf object of POLYGONs with columns `cell_id`, `neighbor_id` and, with
values, `value`, `neighbor_value` and `change`
(`neighbor_value - value`).

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`hex_rays`](https://gillescolling.com/hexify/reference/hex_rays.md),
[`get_neighbors`](https://gillescolling.com/hexify/reference/get_neighbors.md),
[`wall_metrics`](https://gillescolling.com/hexify/reference/wall_metrics.md)

## Examples

``` r
g <- hex_grid(resolution = 3, aperture = 3)
cells <- grid_rect(c(0, 40, 30, 60), g)$cell_id
ctr <- cell_to_lonlat(cells, g)
tri <- hex_triangles(cells, g, value = ctr$lat_deg)
plot(tri["change"], border = NA)
```
