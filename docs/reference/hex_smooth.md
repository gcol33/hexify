# Smooth cell values over neighbouring cells

Replaces each cell's value by the weighted mean of its own value and the
values of its edge neighbours, and repeats this `steps` times: the
neighbour-matrix smoother of Carr et al. (1997). Each neighbour has
weight 1 and the cell itself `self_weight`. On an equal-area grid every
cell has the same area, so the mean is also the area-weighted mean.

## Usage

``` r
hex_smooth(cell_id, value, grid, steps = 1L, self_weight = 1, group = NULL)
```

## Arguments

- cell_id:

  Cell IDs, each listed once.

- value:

  Numeric values, one per cell.

- grid:

  A HexGridInfo or HexData object.

- steps:

  Number of smoothing steps. Each step widens the reach by one ring of
  cells.

- self_weight:

  Weight of a cell's own value relative to each neighbour's. `0`
  replaces a value by the mean of its neighbours.

- group:

  Optional vector, one entry per cell; only neighbours in the same group
  contribute.

## Value

Numeric vector of smoothed values, in the order of `cell_id`.

## Details

Only the cells in `cell_id` take part: a neighbour that is not listed,
or whose value is `NA`, does not contribute, so a cell at the edge of
the data averages over the neighbours it has. A pentagon averages over
its five neighbours. A cell whose value is `NA` stays `NA`.

`group` restricts the smoothing to cells of the same group, such as land
and ocean: a neighbour in another group does not contribute, so values
do not mix across a coastline. A cell whose group is `NA` keeps its
value and contributes to no other cell.

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md)
for the neighbours a step averages over

## Examples

``` r
# \donttest{
g <- hex_grid(resolution = 4, aperture = 3)
cells <- grid_global(g)
set.seed(1)
noise <- rnorm(nrow(cells))
smooth <- hex_smooth(cells$cell_id, noise, g, steps = 3)
c(sd(noise), sd(smooth))

# Keep land and ocean apart
ctr <- sf::st_as_sf(cell_to_lonlat(cells$cell_id, g),
                    coords = c("lon_deg", "lat_deg"), crs = 4326)
land <- lengths(sf::st_intersects(ctr, hexify_world)) > 0
coastal <- hex_smooth(cells$cell_id, noise, g, steps = 3, group = land)
# }
```
