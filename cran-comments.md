## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission (v0.6.4)

Fixes for issues flagged in v0.6.3 incoming checks:

1. Removed compiled object files (`.o`) from source tarball — caused
   installation ERROR on Debian and NOTE on all platforms
2. Wrapped `plot_globe()` examples in `\donttest{}` — was 608s on win-builder
3. Reworded DESCRIPTION to avoid "vendored" spelling flag

## Changes in v0.6.3 (carried forward)

1. `cell_to_sf()` now applies `st_wrap_dateline()` automatically for correct
   antimeridian rendering on flat map projections
2. Fixed antimeridian gaps in `plot_globe()` orthographic projection
3. Fixed `sprintf` → `snprintf` in H3 C code (compiled code WARNING)
4. `as_sf(geometry = "polygon")` routes all grids through `cell_to_sf()` for
   consistent antimeridian handling

## Test environments

* local Windows 11 / WSL2, R 4.5.2
* GitHub Actions (ubuntu-latest release/devel/oldrel, macOS-latest, windows-latest)

## Downstream dependencies

There are currently no downstream dependencies for this package.
