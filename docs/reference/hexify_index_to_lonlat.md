# Convert index string to longitude/latitude

Returns the cell center coordinates for a given index.

## Usage

``` r
hexify_index_to_lonlat(
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

  Index encoding. Default "auto" infers from aperture.

## Value

Named numeric vector with lon and lat in degrees

## Examples

``` r
coords <- hexify_index_to_lonlat("051223", aperture = 3)
```
