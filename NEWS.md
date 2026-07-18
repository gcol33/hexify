# hexify 0.7.4 (development version)

## New features

* Hierarchical navigation now works for mixed aperture `"4/3"` (ISEA43H) grids:
  `get_parent()`, `get_children()`, and `cell_to_index()` no longer error on
  these grids (#31). Mixed 4/3 grids have no DGGRID-standard hierarchical index,
  so hexify defines the hierarchy geometrically -- a cell's parent is the
  coarser cell whose lattice point contains the cell's centre. This is the
  correct relationship for how these grids quantise (a single scaled
  quantisation per resolution, not a step-composed subdivision), and it
  round-trips `cell -> cell_to_index() -> cell` exactly, including the seam and
  icosahedron-vertex cells that a purely local walk would miss.

# hexify 0.7.3

**Bug fixes, package hygiene, and test coverage**

## Fixes

* Fixed `hexify_lonlat_to_index()`/`hexify_index_to_lonlat()` for aperture 3:
  quantization skipped the quad-frame fold that aperture 4/7 already had,
  producing indices inconsistent with `lonlat_to_cell()`/`cell_to_index()`
  for essentially all points.
* Fixed `is_pentagon()` for ISEA aperture 7, which undercounted or threw at
  resolution >= 1; it now decodes each cell's own (i, j) instead of relying
  on an unreliable forward computation.
* Fixed `hex_compact()` emitting duplicate cell IDs when the input already
  contained a parent cell alongside all 7 of its children.
* Fixed `hexify_assign()`'s Z3 backend discarding the sign of quantized
  (i, j), causing sign-flipped points to collide on the same cell ID.
* Fixed NaN/Inf inputs reaching undefined behavior in `snyder_forward()`'s
  sort comparator and silently returning garbage from hex quantization
  instead of erroring.
* `cell_to_index()`/`get_parent()`/`get_children()` no longer produce a
  confusing internal bound error (or, with a naive fix, silently wrong cell
  IDs) for mixed aperture `"4/3"` grids. In 0.7.3 these raised a clear
  "not implemented" error; 0.7.4 implements the navigation geometrically
  (see above).
* `hex_extract()` now respects `cells=`/`boundary=` when `grid` is a
  HexData object (previously silently ignored).
* `hex_browse()`'s data.frame input mode no longer crashes on duplicate
  `cell_id` rows.
* `plot_globe()`'s `resolve_center()` now validates a named `center`
  vector's names instead of silently building an invalid PROJ string.
* Hierarchical-index functions (`hexify_cell_to_index()` and siblings) now
  validate `resolution`/`aperture` like `hexify_lonlat_to_cell()` already
  did; the C++ layer now enforces the max resolution for aperture 3/4 the
  same way it already did for aperture 7.
* `as_dggrid()` now accepts a modern `HexGridInfo` object, matching
  `dgverify()`.
* `HexData`'s validity check no longer has a blind spot for empty
  `cell_id`/`cell_center` paired with non-empty `data`.
* `hex_grid()` and legacy `hexify_grid()` now give clear errors for
  non-numeric/NA/non-positive `resolution`, `crs`, or `area`, instead of
  base-R errors or a silent `resolution = NaN`.
* `plot_globe(exclude_antarctica = TRUE)` now warns instead of silently
  no-op'ing for custom `land_data` without a recognized country-name
  column.
* `import_h3()` now drops NA cell IDs (with a warning) when `data` is
  attached, matching `hexify()`'s existing NA-coordinate handling.
* `prepare_fill_column()` now warns when `breaks=` is supplied for a
  discrete value column instead of silently discarding it.

## Package hygiene

* Removed dead code (`index_to_cell_internal()`), an orphaned test
  fixture, and a committed rendered vignette; fixed a native-pipe usage
  that required a newer R than the package declares; fixed stale roxygen
  text; documented `HexData[`'s `drop = FALSE` default explicitly.
* C++: replaced a raw `malloc`/`free` digit buffer with `std::vector`;
  an out-of-range digit now throws instead of silently clamping; added a
  resolution/string-length consistency check to `z3::decode()`;
  `get_children_indices()` no longer swallows unrelated exceptions behind
  a catch-all around the max-resolution boundary check.

## Tests

* Strengthened ~20 `hexify_heatmap()` tests that previously only checked
  the return type is a ggplot object to assert on the actual plot
  data/mapping (fill column, scale type, colors, bins, limits, labels).
* Added coverage for the `Raster*`/`SpatRaster` basemap path and for
  `hex_browse()`'s value-to-fill-color mapping.

# hexify 0.7.2

* Fixed aperture 4/7 lon/lat-to-index and index-to-lon/lat conversion:
  quantization was operating on raw icosahedron triangle coordinates
  instead of folding into the quad frame first, producing wrong indices
  for any point on triangles numbered 12-19.
* Implemented `cpp_hex_index_z3_quantize_digits()`, `cpp_hex_index_z3_center()`,
  and the Z3 corners helper, which previously returned placeholder zeros
  instead of real quantized digits/coordinates.
* Fixed `hex_zonal()` row misalignment: results are now keyed off the
  deduplicated `hex_sf$cell_id` order instead of the pre-dedup input,
  and `cells` input now drops `NA`/duplicate values before lookup.
* Fixed `grid_clip()` argument order in `hex_zonal()`'s boundary path.

# hexify 0.7.1

**Code quality and documentation**

## Improvements

* Refactored `rcpp_aperture.cpp`: extracted parameterized helpers for
  quantize/center/corners, eliminating copy-paste across 3 apertures.
* Consolidated input validation: conversion functions now use shared
  `validate_aperture()` / `validate_resolution()` from `constants.R`.
* Optimized `is_pentagon()` for ISEA grids: computes pentagon cell IDs
  directly from quad coordinates instead of looping through the lon/lat
  pipeline.
* Stricter aperture parsing in `hexify()`: rejects unexpected values with
  an informative error instead of silently coercing via `as.integer()`.
* Deduplicated `theme_minimal()` calls in plot methods with `.theme_clean()`
  helper.

## Documentation

* Added "Why Hexagonal Grids?" section to README (equal area, uniform
  adjacency, low shape distortion).
* Added "Known Limitations" section to README (H3 resolution cap, pentagons,
  projection precision).

## Tests

* New `test-edge-cases.R` with 32 tests covering poles, antimeridian,
  dateline, equator, roundtrip stability, pentagon invariants, neighbor
  symmetry, and resolution-0 cell counts.

# hexify 0.7.0

**Spatial analysis primitives**

## Hotfix

* Fixed aperture 7 cell encoding to use Class III substrate quantization.
* Updated `max_cell_id()` for aperture 7 bounding box.

## New features

* New `get_neighbors()`: returns k-ring (disk) of neighboring cells for
  both ISEA and H3 grids. Supports `k > 1` for multi-ring expansion,
  `distances = TRUE` for ring distance output, and vectorized cell input.
  ISEA backend uses axial coordinate offsets with lon/lat fallback for
  cross-quad boundaries. H3 backend uses vendored `gridDisk` /
  `gridDiskDistances` / `gridRingUnsafe`.
* New `hex_summarize()`: cell-level data aggregation with tidyeval support.
  Groups by cell_id, applies user-defined summary expressions, returns a
  data.frame with cell centers, areas, and point counts. Supports
  `geometry = TRUE` for sf output.
* New C++ bindings: `cpp_h3_gridDisk`, `cpp_h3_gridDiskDistances`,
  `cpp_h3_gridRingUnsafe`, `cpp_get_neighbors_isea`, `cpp_get_neighbors_z7`

# hexify 0.6.5

* Fixed empty translation unit warning in vendored H3 `h3Assert.c`
  (clang 21 `-Wempty-translation-unit`)

# hexify 0.6.4

* Removed compiled object files (`.o`) from source tarball that caused
  installation failure on Linux (Debian) and NOTE on all platforms
* Wrapped `plot_globe()` examples in `\donttest{}` to reduce check time
  (was 608s on win-builder)
* Reworded DESCRIPTION to avoid "vendored" spelling flag

# hexify 0.6.3

* `cell_to_sf()` now applies `sf::st_wrap_dateline()` automatically, fixing
  horizontal streaks on flat map projections (Plate Carrée, Robinson, etc.) for
  hexagons crossing the ±180° antimeridian
* `as_sf(x, geometry = "polygon")` now routes all grids (ISEA and H3) through
  `cell_to_sf()`, ensuring consistent antimeridian handling
* `hexify_cell_to_sf()` gains antimeridian normalization matching `cell_to_sf()`
* Vignettes no longer require manual `st_wrap_dateline()` calls

# hexify 0.6.2

**Native H3 backend — zero external dependencies**

* Vendored H3 v4.4.1 C source, replacing the `h3o` R package dependency
* H3 is now always available — no optional install, no Suggests
* Uses H3 experimental polygon fill for full spatial coverage in `hexify()`
* New native C++ bindings: `h3_lat_lng_to_cell`, `h3_cell_to_boundary`,
  `h3_cell_to_parent`, `h3_cell_to_children`, `h3_polygon_to_cells`,
  `h3_cell_area_km2`, `h3_cell_to_lat_lng`

# hexify 0.6.1

* Warn when `aperture` is passed with `type = "h3"` (ignored parameter)
* Expanded H3 resolution guidance in `hex_grid()` documentation
* Extended H3 vignette resolution table to full range (0-15)
* Added `h3_crosswalk()` example to H3 vignette

# hexify 0.6.0

**ISEA–H3 crosswalk and per-cell area**

* New `h3_crosswalk()`: bidirectional mapping between ISEA and H3 cell IDs,
  with automatic resolution matching and per-cell area comparison
* New `cell_area()`: returns geodesic area (km²) for each cell — constant for
  ISEA (equal-area), location-dependent for H3, with session-scoped caching
* HexData `$cell_area_km2`, `[["cell_area_km2"]]`, and `as.data.frame()` now
  return per-cell areas for H3 grids instead of the grid-wide average
* Internal: extracted `closest_h3_resolution()` helper shared by `hex_grid()`
  and `h3_crosswalk()`

# hexify 0.5.0

**H3 grid support**

* Added H3 (Uber) as a first-class grid type: `hex_grid(resolution = 8, type = "h3")`
* All core functions work with H3 grids: `hexify()`, `cell_to_sf()`, `grid_rect()`,
  `grid_global()`, `grid_clip()`, `get_parent()`, `get_children()`
* H3 support requires the `h3o` package (Suggests, not required for ISEA workflows)
* New `hexify_compare_resolutions(type = "h3")` for H3 resolution table
* `dgearthstat()` now accepts HexGridInfo objects directly
* New `grid_type` slot on HexGridInfo: `"isea"` (default) or `"h3"`
* HexData `cell_id` slot supports character (H3) and numeric (ISEA) cell IDs
* Backward compatible: all existing ISEA workflows unchanged

# hexify 0.3.10

**Hotfix for geometry issues**

* Fixed invalid pentagon geometries that caused gaps in global grids
* Fixed antimeridian-crossing polygons using `st_wrap_dateline()`
* Added polar cap sampling to `grid_global()` to include cells above ±85° latitude

# hexify 0.3.6

* Reduced test suite runtime for CRAN by skipping detailed consistency tests
  (full tests still run locally via NOT_CRAN=true)
* Fixed CRAN incoming check NOTE: "Overall checktime 15 min > 10 min"

# hexify 0.3.5

* Simplified plot examples to reduce runtime
* Added non-standard files to .Rbuildignore
* Fixed slow example NOTE

# hexify 0.3.4

**Hotfix for CRAN UBSAN check failure**

* Fixed undefined behavior in Snyder ISEA projection causing NaN values
  (UBSAN error on M1 Mac CRAN check: "nan is outside the range of representable
  values of type 'long long'" at coordinate_transforms.cpp:243-244)
* Root cause: floating-point precision in face assignment could project points
  onto geometrically invalid triangle faces near icosahedron edges
* Solution: added projection validation with face-retry logic - if validation
  fails (z > DH), automatically tries adjacent faces until valid

# hexify 0.3.3

* CRAN submission fixes

# hexify 0.3.2

* Documentation improvements
* Fixed hierarchical index functions
* Increased test coverage to 90%

# hexify 0.3.1

* Minor bug fixes

# hexify 0.3.0

* Major refactoring of coordinate transformation system
* Improved performance for large grids

# hexify 0.2.0

* Added dggridR compatibility functions (`as_dggrid()`, `from_dggrid()`)
* New hierarchical indexing support with H-index functions
* Improved coordinate conversion pipeline
* Added grid statistics functions

# hexify 0.1.0

* Initial release
* ISEA discrete global grid implementation
* Support for apertures 3, 4, 7, and mixed 4/3
* Compatible with dggridR output
