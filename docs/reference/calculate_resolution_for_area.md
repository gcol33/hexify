# Calculate resolution for target area

Uses the cell count formula N = d \* aperture^res + 2, d the solid's
diamond quads; on the icosahedron (d = 10) this is the
'ISEA3H'/'ISEA4H'/'ISEA7H' count, which matches 'dggridR' resolution
numbering exactly.

## Usage

``` r
calculate_resolution_for_area(
  target_area_km2,
  aperture = 3,
  radius_km = EARTH_RADIUS_KM,
  polyhedron = "icosahedron"
)
```

## Arguments

- target_area_km2:

  Target area in square kilometers

- aperture:

  Aperture (3, 4, or 7)

- radius_km:

  Radius of the body, in kilometers

- polyhedron:

  The solid the grid is built on

## Value

Resolution level, not rounded. A target larger than the cells of
resolution 0 gives -Inf.
