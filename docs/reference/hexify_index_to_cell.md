# Convert index string to cell coordinates

Decodes a hierarchical index string back to cell coordinates.

## Usage

``` r
hexify_index_to_cell(
  index,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- index:

  Index string

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding used. Default "auto" infers from aperture.

## Value

List with face, i, j, and resolution

## Examples

``` r
cell <- hexify_index_to_cell("051223", aperture = 3)
```
