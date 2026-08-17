# Get Neighboring Cells

Returns the k-ring (disk) of cells neighboring the input cells. For
`k = 1`, returns the immediate 6 neighbors (5 for pentagons). For
`k > 1`, returns all cells within `k` grid hops.

## Usage

``` r
get_neighbors(cell_id, grid, k = 1L, include_self = FALSE, distances = FALSE)
```

## Arguments

- cell_id:

  Cell IDs to find neighbors for. Numeric vector for ISEA grids,
  character vector for H3 grids.

- grid:

  A HexGridInfo or HexData object specifying the grid.

- k:

  Integer. Ring distance (default 1). `k = 1` returns immediate
  neighbors, `k = 2` includes neighbors-of-neighbors, etc.

- include_self:

  Logical. If `TRUE`, include the input cell in the result (default
  `FALSE`).

- distances:

  Logical. If `TRUE`, return a data.frame with cell IDs and their ring
  distance from the origin (default `FALSE`).

## Value

If `distances = FALSE` (default): a list of cell ID vectors, one per
input cell. If `distances = TRUE`: a list of data.frames with columns
`cell_id` and `ring_distance`.

## Details

For **ISEA grids**, neighbors are computed using axial coordinate
offsets in the quad IJ space. Cells at quad boundaries are handled by
reprojection through lon/lat coordinates.

For **H3 grids**, neighbors use the vendored H3 `gridDisk` /
`gridDiskDistances` functions.

Pentagon cells (the 12 icosahedron vertices) have only 5 neighbors
instead of the usual 6.

## See also

[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md) for
creating HexData objects,
[`hex_distance()`](https://gillescolling.com/hexify/reference/hex_distance.md)
for grid distances between cells

## Examples

``` r
# \donttest{
# ISEA grid neighbors
g <- hex_grid(area_km2 = 1000)
cell <- lonlat_to_cell(10, 50, g)
nbrs <- get_neighbors(cell, g)
nbrs[[1]]

# H3 grid neighbors
g_h3 <- hex_grid(resolution = 5, type = "h3")
cell_h3 <- lonlat_to_cell(10, 50, g_h3)
get_neighbors(cell_h3, g_h3, k = 2)

# With distances
get_neighbors(cell, g, k = 2, distances = TRUE)
# }
```
