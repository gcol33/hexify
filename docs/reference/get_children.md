# Get children cells

Returns the child cells at a finer resolution.

## Usage

``` r
get_children(cell_id, grid, levels = 1L, as_sf = FALSE)
```

## Arguments

- cell_id:

  Cell IDs: integer64 for ISEA grids (whole numbers below 2^53 and digit
  strings are accepted), character for H3 grids

- grid:

  A HexGridInfo or HexData object

- levels:

  Number of levels down (default 1)

- as_sf:

  If `TRUE`, return the children as sf polygons, one row per child with
  its `parent_id`.

## Value

List of vectors of child cell IDs, one per input cell; with
`as_sf = TRUE`, an sf object with columns `parent_id`, `cell_id` and
`geometry`.

## Examples

``` r
grid <- hex_grid(resolution = 3, aperture = 4)
parent <- lonlat_to_cell(10, 50, grid)
get_children(parent, grid)
kids <- get_children(parent, grid, levels = 2, as_sf = TRUE)
plot(sf::st_geometry(kids))
plot(sf::st_geometry(cell_to_sf(parent, grid)), border = "red", add = TRUE)
```
