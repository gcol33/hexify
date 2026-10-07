# Convert cell coordinates to index string

Converts cell coordinates (`face`, `i`, `j`) to a hierarchical index
string. The index type is selected from `aperture` when
`index_type = "auto"`.

## Usage

``` r
hexify_cell_to_index(
  face,
  i,
  j,
  resolution,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- face:

  Face/quad numbers (0-19)

- i:

  I coordinates

- j:

  J coordinates

- resolution:

  Resolution level

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding: "auto" (default), "z3", "z7", or "zorder"

## Value

A character vector of index strings, one per cell. `face`, `i` and `j`
are read in step and must be the same length.

## Details

Default index types by aperture:

- Aperture 3: Z3 (optimized digit selection)

- Aperture 4: Z-order (Morton curve)

- Aperture 7: Z7 (one base-7 child digit per resolution)

A Z7 index has the form `BBd1...dr`, where `BB` is the two-digit base
cell (00–11), `r` is the resolution, and every child digit is in 0–6. It
is IGEO7's Z7 index, the string DGGRID writes for the same cell, and
decoding and re-encoding a valid index returns the same string. `BB` is
the pentagon the cell descends from, which need not be `face`. Under
each pentagon one direction is missing: a string whose first nonzero
digit is 2 (base cells 00–05) or 5 (06–11) names no cell.

## See also

Other hierarchical index:
[`hexify_compare_indices()`](https://gillescolling.com/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gillescolling.com/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)

## Examples

``` r
idx <- hexify_cell_to_index(5, 10, 15, resolution = 3, aperture = 3)

# Aperture-7 indices contain a two-digit base cell and one digit per level
z7 <- hexify_cell_to_index(5, 1, 1, resolution = 2, aperture = 7)
nchar(z7) == 4L
```
