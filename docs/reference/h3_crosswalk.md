# Crosswalk Between ISEA and H3 Cell IDs

Maps cell IDs between ISEA (equal-area) and H3 grid systems by looking
up each cell's center coordinate in the target grid. This enables
workflows where analysis is done in ISEA (exact equal-area) and
reporting in H3 (industry-standard).

## Usage

``` r
h3_crosswalk(
  cell_id = NULL,
  grid,
  h3_resolution = NULL,
  isea_grid = NULL,
  direction = c("isea_to_h3", "h3_to_isea")
)
```

## Arguments

- cell_id:

  Cell IDs to translate. Numeric for ISEA, character for H3. When `grid`
  is a HexData object and `cell_id` is `NULL`, all cell IDs from the
  data are used.

- grid:

  A HexGridInfo or HexData object. For `direction = "isea_to_h3"`, this
  must be an ISEA grid. For `direction = "h3_to_isea"`, this must be an
  H3 grid.

- h3_resolution:

  Target H3 resolution for `"isea_to_h3"`, or the source H3 resolution
  for `"h3_to_isea"`. When `NULL` (default), the closest H3 resolution
  matching the ISEA cell area is selected automatically.

- isea_grid:

  A HexGridInfo for the target ISEA grid. Required when
  `direction = "h3_to_isea"`.

- direction:

  One of `"isea_to_h3"` (default) or `"h3_to_isea"`.

## Value

A data frame with columns:

- isea_cell_id:

  ISEA cell ID (numeric)

- h3_cell_id:

  H3 cell ID (character)

- isea_area_km2:

  Area of the ISEA cell in km2

- h3_area_km2:

  Geodesic area of the H3 cell in km2

- area_ratio:

  Ratio of ISEA area to H3 area

## Details

The crosswalk works by computing the center coordinate of each source
cell, then finding which cell in the target grid contains that center.
This is a many-to-one mapping: multiple ISEA cells may map to the same
H3 cell (or vice versa) depending on the relative resolutions.

When `h3_resolution` is `NULL` and `direction = "isea_to_h3"`, the H3
resolution whose average cell area is closest to the ISEA cell area is
chosen automatically. This gives the best 1:1 correspondence.

## See also

[`cell_area`](https://gcol33.github.io/hexify/reference/cell_area.md)
for per-cell area computation,
[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
creating grids

## Examples

``` r
# \donttest{
# ISEA -> H3
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
xwalk <- h3_crosswalk(cells, grid)
head(xwalk)

# H3 -> ISEA
h3 <- hex_grid(resolution = 5, type = "h3")
h3_cells <- lonlat_to_cell(c(0, 10), c(45, 50), h3)
xwalk2 <- h3_crosswalk(h3_cells, h3, isea_grid = grid, direction = "h3_to_isea")
# }
```
