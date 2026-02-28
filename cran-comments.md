## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes in v0.6.3

1. `cell_to_sf()` now applies `st_wrap_dateline()` automatically for correct
   antimeridian rendering on flat map projections
2. Fixed antimeridian gaps in `plot_globe()` orthographic projection
3. Fixed `sprintf` → `snprintf` in vendored H3 C code (compiled code WARNING)
4. `as_sf(geometry = "polygon")` routes all grids through `cell_to_sf()` for
   consistent antimeridian handling

## Test environments

* local Windows 11 / WSL2, R 4.5.2
* GitHub Actions (ubuntu-latest release/devel/oldrel, macOS-latest, windows-latest)

## Downstream dependencies

There are currently no downstream dependencies for this package.
