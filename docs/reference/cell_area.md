# Compute per-cell area in km²

Returns the area of each cell in square kilometers. For ISEA grids, all
cells have the same area (equal-area property). For H3 grids, each cell
has a different geodesic area depending on its location.

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

For ISEA grids the area is constant across all cells and is read
directly from the grid specification.

For H3 grids the area varies by latitude. This function computes
geodesic area via
[`sf::st_area()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
on H3 cell polygons, with results cached in a session-scoped environment
so repeated calls for the same cells are fast.

## See also

[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
grid specifications,
[`h3_crosswalk`](https://gcol33.github.io/hexify/reference/h3_crosswalk.md)
for ISEA/H3 interoperability

## Examples

``` r
# ISEA: constant area
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
