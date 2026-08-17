# Import External H3 Cell IDs into hexify

Ingests H3 cell IDs from an external source (another H3 library, a
database, or a CSV file) into hexify. Validates cell IDs, infers the H3
resolution, and optionally attaches data to build a HexData object.

## Usage

``` r
import_h3(cell_ids, data = NULL, validate = TRUE, radius_km = EARTH_RADIUS_KM)
```

## Arguments

- cell_ids:

  Character vector of H3 cell ID strings

- data:

  Optional data frame to attach. Must have the same number of rows as
  `cell_ids`. If `NULL`, returns a HexGridInfo object.

- validate:

  If `TRUE` (default), checks that all cell IDs are valid H3 cells at
  the same resolution before proceeding. Set to `FALSE` to skip
  validation when cell IDs are known to be correct.

- radius_km:

  Radius of the body the cells cover, in kilometers, or a body name
  (default Earth). See
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md).

## Value

If `data = NULL`, a HexGridInfo object for the inferred H3 resolution.
If `data` is provided, a HexData object with data attached at the
specified cells.

## Details

For converting between grid specs, use
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)`(type = "h3")`
directly. For cell-level ISEA/H3 mapping, use
[`h3_crosswalk`](https://gillescolling.com/hexify/reference/h3_crosswalk.md).

H3 cell IDs encode their resolution in the index itself, so no
resolution argument is needed. The resolution is inferred automatically.
All cell IDs must share the same resolution; mixed resolutions produce
an error.

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
creating grids directly,
[`h3_crosswalk`](https://gillescolling.com/hexify/reference/h3_crosswalk.md)
for cell-level ISEA/H3 mapping

## Examples

``` r
# \donttest{
# Import external H3 cell IDs (grid spec only)
h3_ids <- c("8528342bfffffff", "85283473fffffff", "85283447fffffff")
grid <- import_h3(h3_ids)
grid

# Import with data attached
df <- data.frame(species = c("oak", "pine", "birch"), count = c(10, 5, 3))
hd <- import_h3(h3_ids, data = df)
hd
# }
```
