# Assign hex cells ('ISEA3H', aperture 3) for lon/lat

Forward -\> quantize (Z3) -\> center (face) -\> inverse to lon/lat.
Optionally return polygons (sf), in which case sf must be installed.

## Usage

``` r
hexify_assign(
  lon,
  lat,
  effective_res,
  match_dggrid_parity = TRUE,
  make_polygons = FALSE
)
```

## Arguments

- lon, lat:

  numeric vectors (same length), degrees.

- effective_res:

  integer effective resolution (\>= 1).

- match_dggrid_parity:

  logical; TRUE matches 'ISEA3H' parity used by 'dggridR'.

- make_polygons:

  logical; if TRUE, return an sf with hex polygons.

## Value

data.frame with id, face, effective_res, center_lon, center_lat; if
make_polygons=TRUE, an sf with geometry column.
