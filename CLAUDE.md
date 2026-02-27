# hexify Development Notes

## Grid Coverage

hexify supports **all** major hexagonal DGGS through two backends:

- **ISEA** (built-in C++): apertures 3, 4, 7, and mixed 4/3 — resolutions 0-30
- **H3** (via h3o): fixed aperture 7 — resolutions 0-15

This covers every hexagonal grid system that matters:
- ISEA3H, ISEA4H, ISEA7H, ISEA43H = ISEA backend with different aperture settings
- H3 = H3 backend
- IGEO7/Z7 (2025, Sahr) = equal-area aperture-7 hex grid with Z7 indexing = already covered by `hex_grid(aperture = 7)` (ISEA7H with Z7 index). hexify uses the same Z7 hierarchical indexing (7-digit encoding, 0-6 per level) in `src/index_z7.cpp`.
- OpenEAGGR ISEA3H = already covered by `hex_grid(aperture = 3)`
- rHEALPix = diamond-based, not hexagonal — out of scope

No additional grid backends needed.

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
