# Compute per-cell area in km²

Returns the area of each cell in square kilometers. On the ISEA
projection every hexagon of a resolution has the same area, the 12
pentagons of the icosahedron 5/6 of it and the 6 squares of the
octahedron 4/6 of it. On the Fuller projection and for H3 grids, the
area varies from cell to cell.

## Usage

``` r
cell_area(cell_id = NULL, grid)
```

## Arguments

- cell_id:

  Cell IDs to compute area for. For ISEA grids, these are numeric; for
  H3 grids, character strings. When `grid` is a HexData object and
  `cell_id` is `NULL`, all cell IDs from the data are used.

- grid:

  A HexGridInfo or HexData object.

## Value

Named numeric vector of areas in km², one per `cell_id`.

## Details

On the ISEA projection, Snyder's equal-area projection gives every
hexagon of a resolution the same area. A vertex cell, centred on a
vertex of the solid where \\k\\ faces meet, covers \\k\\ of a hexagon's
six sixths, so a grid of \\N\\ cells on a body of area \\S\\ has
hexagons of area \\S / (N - 2)\\, pentagons of \\(5/6) S / (N - 2)\\ on
the icosahedron and squares of \\(4/6) S / (N - 2)\\ on the octahedron.
The grid's `area_km2` is the mean, \\S / N\\.

On the Fuller projection (`hex_grid(projection = "fuller")`) cells are
not equal-area. Each cell's area is its solid angle, summed over its
true boundary, times the body's area over \\4\pi\\, so the cells of a
grid add up to the body's area.

For H3 grids the area varies from cell to cell, by about a factor of 2
across the hexagons of a resolution. The vendored 'H3' library computes
each cell's spherical polygon area as a solid angle, which this function
reads on the grid's body, so a grid built with `radius_km` reports that
body's areas.

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
grid specifications,
[`h3_crosswalk`](https://gillescolling.com/hexify/reference/h3_crosswalk.md)
for ISEA/H3 interoperability

## Examples

``` r
# ISEA: one area for every hexagon
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
cell_area(cells, grid)

# H3: area varies by location
# \donttest{
h3 <- hex_grid(resolution = 5, type = "h3")
h3_cells <- lonlat_to_cell(c(0, 0), c(0, 80), h3)
cell_area(h3_cells, h3)  # equator vs polar — different areas
# }
```
