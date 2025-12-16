# Get index resolution

Returns the resolution level encoded in an index string.

## Usage

``` r
hexify_get_resolution(
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

  Index encoding. Default "auto".

## Value

Integer resolution level
