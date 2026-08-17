# Assign hex cells ('ISEA3H', aperture 3) for lon/lat

Forward projection -\> quantize -\> cell centre -\> inverse to lon/lat,
on the aperture-3 pipeline used by
[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md):
points fold into the non-negative quad frame
([`lonlat_to_cell()`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md)),
and the centre comes back from the cell itself
([`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)),
so the returned centre is the centre of the assigned cell.

## Usage

``` r
hexify_assign(lon, lat, effective_res, make_polygons = FALSE)
```

## Arguments

- lon, lat:

  numeric vectors (same length), degrees.

- effective_res:

  integer effective resolution (\>= 1).

- make_polygons:

  logical; if TRUE, return an sf with hex polygons.

## Value

data.frame with id (Z3 index string), face (quad, 0-11), effective_res,
center_lon, center_lat; if make_polygons=TRUE, an sf with geometry
column.
