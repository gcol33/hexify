# Changelog

## hexify 0.6.2

**Native H3 backend — zero external dependencies**

- Vendored H3 v4.4.1 C source, replacing the `h3o` R package dependency
- H3 is now always available — no optional install, no Suggests
- Uses H3 experimental polygon fill for full spatial coverage in
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)
- New native C++ bindings: `h3_lat_lng_to_cell`, `h3_cell_to_boundary`,
  `h3_cell_to_parent`, `h3_cell_to_children`, `h3_polygon_to_cells`,
  `h3_cell_area_km2`, `h3_cell_to_lat_lng`

## hexify 0.6.1

- Warn when `aperture` is passed with `type = "h3"` (ignored parameter)
- Expanded H3 resolution guidance in
  [`hex_grid()`](https://gcol33.github.io/hexify/reference/hex_grid.md)
  documentation
- Extended H3 vignette resolution table to full range (0-15)
- Added
  [`h3_crosswalk()`](https://gcol33.github.io/hexify/reference/h3_crosswalk.md)
  example to H3 vignette

## hexify 0.6.0

**ISEA–H3 crosswalk and per-cell area**

- New
  [`h3_crosswalk()`](https://gcol33.github.io/hexify/reference/h3_crosswalk.md):
  bidirectional mapping between ISEA and H3 cell IDs, with automatic
  resolution matching and per-cell area comparison
- New
  [`cell_area()`](https://gcol33.github.io/hexify/reference/cell_area.md):
  returns geodesic area (km²) for each cell — constant for ISEA
  (equal-area), location-dependent for H3, with session-scoped caching
- HexData `$cell_area_km2`, `[["cell_area_km2"]]`, and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) now
  return per-cell areas for H3 grids instead of the grid-wide average
- Internal: extracted `closest_h3_resolution()` helper shared by
  [`hex_grid()`](https://gcol33.github.io/hexify/reference/hex_grid.md)
  and
  [`h3_crosswalk()`](https://gcol33.github.io/hexify/reference/h3_crosswalk.md)

## hexify 0.5.0

**H3 grid support**

- Added H3 (Uber) as a first-class grid type:
  `hex_grid(resolution = 8, type = "h3")`
- All core functions work with H3 grids:
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md),
  [`cell_to_sf()`](https://gcol33.github.io/hexify/reference/cell_to_sf.md),
  [`grid_rect()`](https://gcol33.github.io/hexify/reference/grid_rect.md),
  [`grid_global()`](https://gcol33.github.io/hexify/reference/grid_global.md),
  [`grid_clip()`](https://gcol33.github.io/hexify/reference/grid_clip.md),
  [`get_parent()`](https://gcol33.github.io/hexify/reference/get_parent.md),
  [`get_children()`](https://gcol33.github.io/hexify/reference/get_children.md)
- H3 support requires the `h3o` package (Suggests, not required for ISEA
  workflows)
- New `hexify_compare_resolutions(type = "h3")` for H3 resolution table
- [`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md)
  now accepts HexGridInfo objects directly
- New `grid_type` slot on HexGridInfo: `"isea"` (default) or `"h3"`
- HexData `cell_id` slot supports character (H3) and numeric (ISEA) cell
  IDs
- Backward compatible: all existing ISEA workflows unchanged

## hexify 0.3.10

CRAN release: 2026-02-04

**Hotfix for geometry issues**

- Fixed invalid pentagon geometries that caused gaps in global grids
- Fixed antimeridian-crossing polygons using
  [`st_wrap_dateline()`](https://r-spatial.github.io/sf/reference/st_transform.html)
- Added polar cap sampling to
  [`grid_global()`](https://gcol33.github.io/hexify/reference/grid_global.md)
  to include cells above ±85° latitude

## hexify 0.3.6

- Reduced test suite runtime for CRAN by skipping detailed consistency
  tests (full tests still run locally via NOT_CRAN=true)
- Fixed CRAN incoming check NOTE: “Overall checktime 15 min \> 10 min”

## hexify 0.3.5

- Simplified plot examples to reduce runtime
- Added non-standard files to .Rbuildignore
- Fixed slow example NOTE

## hexify 0.3.4

**Hotfix for CRAN UBSAN check failure**

- Fixed undefined behavior in Snyder ISEA projection causing NaN values
  (UBSAN error on M1 Mac CRAN check: “nan is outside the range of
  representable values of type ‘long long’” at
  coordinate_transforms.cpp:243-244)
- Root cause: floating-point precision in face assignment could project
  points onto geometrically invalid triangle faces near icosahedron
  edges
- Solution: added projection validation with face-retry logic - if
  validation fails (z \> DH), automatically tries adjacent faces until
  valid

## hexify 0.3.3

CRAN release: 2026-01-20

- CRAN submission fixes

## hexify 0.3.2

- Documentation improvements
- Fixed hierarchical index functions
- Increased test coverage to 90%

## hexify 0.3.1

- Minor bug fixes

## hexify 0.3.0

- Major refactoring of coordinate transformation system
- Improved performance for large grids

## hexify 0.2.0

- Added dggridR compatibility functions
  ([`as_dggrid()`](https://gcol33.github.io/hexify/reference/as_dggrid.md),
  [`from_dggrid()`](https://gcol33.github.io/hexify/reference/from_dggrid.md))
- New hierarchical indexing support with H-index functions
- Improved coordinate conversion pipeline
- Added grid statistics functions

## hexify 0.1.0

- Initial release
- ISEA discrete global grid implementation
- Support for apertures 3, 4, 7, and mixed 4/3
- Compatible with dggridR output
