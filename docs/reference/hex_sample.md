# Random points inside cells

Draws points uniformly over the area of each cell on the sphere, as
DGGRID's `randpts` output does. Each point is assigned back to its cell
by
[`lonlat_to_cell`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md),
so the points of a cell are exactly the points the grid places in it.
Use [`set.seed()`](https://rdrr.io/r/base/Random.html) for a
reproducible draw.

## Usage

``` r
hex_sample(cell_id, grid, n = 1L)
```

## Arguments

- cell_id:

  Cell IDs. integer64 for ISEA grids, character for H3 grids.

- grid:

  A HexGridInfo or HexData object.

- n:

  Number of points per cell: one number, or one per cell.

## Value

A data frame with columns `cell_id`, `lon` and `lat`, the points of each
cell together in the order of `cell_id`.

## Details

Points are drawn uniformly in a spherical cap around the cell centre
that holds the whole cell, and those outside the cell are drawn again.

## See also

[`cell_to_lonlat`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
for cell centres

## Examples

``` r
grid <- hex_grid(resolution = 3, aperture = 3)
cell <- lonlat_to_cell(10, 50, grid)
set.seed(1)
pts <- hex_sample(cell, grid, n = 200)
plot(sf::st_geometry(cell_to_sf(cell, grid)))
points(pts$lon, pts$lat, pch = 20, cex = 0.5)
```
