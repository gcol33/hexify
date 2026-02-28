## R CMD check results

0 errors | 0 warnings | 1 note

The installed size NOTE (~6 MB) is expected: the package includes a C++
core for the ISEA projection engine and the H3 v4.4.1 C library compiled
from source. Debug symbols are now stripped in Makevars to minimize size.

## Resubmission (v0.6.5)

Fixes for issues flagged in v0.6.4 incoming checks:

5. Fixed empty translation unit WARNING in `h3Assert.c` (clang 21
   `-Wempty-translation-unit` on Debian)

Fixes carried from v0.6.4:

1. Removed compiled object files (`.o`) from source tarball — caused
   installation ERROR on Debian and NOTE on all platforms
2. Wrapped `plot_globe()` examples in `\donttest{}` — was 608s on win-builder
3. Reworded DESCRIPTION to avoid "vendored" spelling flag
4. Added `strip -S` in Makevars to reduce installed library size (NOTE)

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
