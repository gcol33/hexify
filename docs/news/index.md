# Changelog

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
