## R CMD check results

0 errors | 0 warnings | 1 note

* This release adds ISEA–H3 crosswalk functions and per-cell area computation.

The NOTE is "Days since last update" which is expected.

## Changes in v0.6.0

1. New `h3_crosswalk()` for bidirectional ISEA ↔ H3 cell ID mapping
2. New `cell_area()` returning geodesic area per cell (constant for ISEA, varies for H3)
3. HexData virtual columns (`$cell_area_km2`) now return per-cell areas for H3 grids

## Test environments

* local Windows 11 / WSL2, R 4.5.2
* GitHub Actions (ubuntu, macOS, windows)

## Downstream dependencies

There are currently no downstream dependencies for this package.
