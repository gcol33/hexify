# Convert Cell ID to PLANE coordinates

Converts DGGRID-compatible cell IDs directly to PLANE coordinates.
Returns the cell center in the unfolded icosahedron layout.

## Usage

``` r
hexify_cell_to_plane(cell_id, resolution, aperture = 3L)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs (1-based)

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Data frame with columns:

- plane_x:

  X coordinate in PLANE space (range ~0 to 5.5)

- plane_y:

  Y coordinate in PLANE space (range ~0 to 1.73)

## Details

Compatible with dggridR's dgSEQNUM_to_PLANE().

## See also

[`hexify_icosa_tri_to_plane`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_plane.md)
for triangle conversion,
[`hexify_lonlat_to_plane`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_plane.md)
for lon/lat conversion

## Examples

``` r
if (FALSE) { # \dontrun{
# Get PLANE coordinates for cells
plane <- hexify_cell_to_plane(cell_id = c(100, 200, 300),
                               resolution = 5, aperture = 3)
plot(plane$plane_x, plane$plane_y)
} # }
```
