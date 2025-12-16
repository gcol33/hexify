# Convert cell coordinates to index string

Converts DGGRID cell coordinates (face, i, j) to a hierarchical index
string. The index type is automatically selected based on aperture
unless specified.

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

  Face/quad number (0-19)

- i:

  I coordinate

- j:

  J coordinate

- resolution:

  Resolution level

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding: "auto" (default), "z3", "z7", or "zorder"

## Value

Index string (e.g., "051223")

## Details

Default index types by aperture:

- Aperture 3: Z3 (optimized digit selection)

- Aperture 4: Z-order (Morton curve)

- Aperture 7: Z7 (hierarchical with Class III handling)

## Examples

``` r
idx <- hexify_cell_to_index(5, 10, 15, resolution = 3, aperture = 3)
```
