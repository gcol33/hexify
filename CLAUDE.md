# hexify Development Notes

## Grid Coverage

hexify supports **all** major hexagonal DGGS through two backends:

- **ISEA** (built-in C++): apertures 3, 4, 7, and any mixed sequence of them — resolutions 0-30.
  `hex_grid(aperture = "4/3")` / `"4/7"` / `"7/4"` name a family (first `floor(res/2)` levels
  take the first aperture), and `aperture = c(4, 4, 7, 3)` names one aperture per level.
- **H3** (via h3o): fixed aperture 7 — resolutions 0-15

This covers every hexagonal grid system that matters:
- ISEA3H, ISEA4H, ISEA7H, ISEA43H = ISEA backend with different aperture settings
- H3 = H3 backend
- IGEO7/Z7 (2025, Sahr) = equal-area aperture-7 hex grid with Z7 indexing = already covered by `hex_grid(aperture = 7)` (ISEA7H with Z7 index). hexify uses the same Z7 hierarchical indexing (7-digit encoding, 0-6 per level) in `src/index_z7.cpp`.
- OpenEAGGR ISEA3H = already covered by `hex_grid(aperture = 3)`
- rHEALPix = diamond-based, not hexagonal — out of scope

No additional grid backends needed.

## Bodies

An ISEA grid is sized on any sphere: `hex_grid(radius_km = )` takes a radius in km
or a body name (`"mars"`, `"moon"`, `"titan"`, ...). The grid object carries the
radius in its `radius_km` slot, and every function reporting kilometres reads it
through `grid_radius_km()` rather than the Earth constant. Cell geometry is
angular and radius-free, so only areas, diagonals, spacings and the
resolution-for-area inversion change. H3 is fixed to Earth's radius by its
vendored C library.

## Vendored H3 C Library

The H3 backend uses vendored C source from Uber's H3 library in `src/h3/`. This is a direct copy of upstream code — **do not rewrite or restyle it**. Benefits:
- Zero external dependencies at install time
- Easy to update: drop in new upstream C files when Uber releases a new version
- Apache 2.0 license (compatible with MIT), kept at `src/h3/LICENSE`

Only make minimal targeted fixes when R CMD check flags specific symbols (e.g., `sprintf` → `snprintf`). Never refactor vendored code for style.

## Build & Check

Use Windows R (not WSL R):
```bash
"/mnt/c/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::document()'
"/mnt/c/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::check(args = "--no-manual")'
```

pkgdown site:
```bash
"/mnt/c/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'source("~/.R/build_pkgdown.R"); build_pkgdown_site()'
```

## Git Push

Use Windows git/gh for remote operations (WSL SSH agent doesn't persist):
```bash
cmd.exe /C "cd /d C:\Users\Gilles Colling\documents\dev\hexify && git push origin main"
cmd.exe /C "cd /d C:\Users\Gilles Colling\documents\dev\hexify && gh release create v0.x.x --title \"title\" --notes \"notes\""
```
