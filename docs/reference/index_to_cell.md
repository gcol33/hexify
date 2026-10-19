# Convert Hierarchical Index to Cell ID

The inverse of
[`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md):
the cells of `grid` that index strings, IGEO7 integer forms or OGC zone
identifiers name.

## Usage

``` r
index_to_cell(
  index,
  grid,
  form = c("string", "int", "hex", "monotonic", "textZIRS", "uint64ZIRS")
)
```

## Arguments

- index:

  Identifiers in the form `form` names: character for `"string"`,
  `"hex"` and `"textZIRS"`, integer64 (or whole numbers below 2^53, or
  for `"uint64ZIRS"` digit strings) for `"int"`, `"monotonic"` and
  `"uint64ZIRS"`

- grid:

  A HexGridInfo or HexData object

- form:

  The form of `index`, as in
  [`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)

## Value

Cell IDs of `grid`: integer64 for ISEA grids, character for H3 grids.
Every index must name a cell of the grid's resolution; an OGC identifier
that names no zone of the grid reads as `NA`.

## Details

A packed or hexadecimal IGEO7 index carries its resolution in its
trailing 7 digits, and an OGC identifier in its level; a monotonic ID is
read at the grid's resolution. See
[`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)
for the forms.

## See also

[`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)

## Examples

``` r
g <- hex_grid(resolution = 4, aperture = 7)
index_to_cell("033000", g)
index_to_cell("3600ffffffffffff", g, form = "hex")
index_to_cell(bit64::as.integer64("8232"), g, form = "monotonic")

ogc <- hex_grid(resolution = 8, aperture = 3, ellipsoid = "WGS84",
                orientation = "ogc")
cell_to_lonlat(index_to_cell("E6-317-A", ogc, form = "textZIRS"), ogc)
```
