# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

hexify is a minimal R package for assigning geographic points to
equal-area hexagonal cells using the ISEA discrete global grid system.
Produces output identical to dggridR with a simpler API. Supports
apertures 3, 4, 7, and mixed 4/3.

**Primary function:** `hexify(df, lon, lat, area)` → returns input with
`cell_id`, `cell_cen_lon`, `cell_cen_lat`, `cell_area`, `cell_diag`
columns.

## Build & Development Commands

Use `Rscript` for command-line execution (not `R -e`).

``` bash
Rscript -e "devtools::load_all()"                              # Load for dev
Rscript -e "devtools::test()"                                  # Run all tests
Rscript -e "testthat::test_file('tests/testthat/test-hexify.R')" # Single test
Rscript -e "devtools::document()"                              # Rebuild docs
Rscript -e "Rcpp::compileAttributes()"                         # Rebuild C++ bindings
```

## Architecture

### Coordinate Pipeline

    lon/lat → Snyder forward → (face, tx, ty)   [Icosa Triangle]
           → icosa_tri_to_quad_ij → (quad, i, j) [Quad IJ]
           → quad_ij_to_cell → SEQNUM            [cell ID]

### C++ Layer (src/)

| File | Purpose |
|----|----|
| `icosahedron.cpp` | Face geometry, `which_face()` |
| `projection_forward.cpp` | lon/lat → (face, tx, ty) |
| `projection_inverse.cpp` | (face, x, y) → lon/lat (Newton-Raphson) |
| `aperture.cpp` | Quantization for apertures 3, 4, 7 |
| `aperture_sequence.cpp` | Mixed aperture (4/3) support |
| `cell_numbering.cpp` | Icosa Triangle ↔︎ Quad IJ ↔︎ cell ID conversions |
| `coordinate_transforms.cpp` | Coordinate space utilities |
| `rcpp_*.cpp` | R bindings |

### R Layer (R/)

| File                  | Purpose                              |
|-----------------------|--------------------------------------|
| `hexify.R`            | Main entry point                     |
| `hexify_grid.R`       | Grid construction, resolution lookup |
| `hexify_projection.R` | Projection wrappers                  |
| `hexify_polygons.R`   | Polygon generation                   |

## Coordinate Systems

| System | Components | Description |
|----|----|----|
| GEO | lon, lat | WGS84 degrees |
| Icosa Triangle | icosa_triangle_face (0-19), icosa_triangle_x, icosa_triangle_y | Triangle face + projected coords |
| Quad XY | quad (0-11), quad_x, quad_y | Quad + continuous coordinates |
| Quad IJ | quad (0-11), i, j | Quad + integer grid indices |
| SEQNUM/Cell ID | single integer | Global cell ID (dggridR compatible) |

Triangle-to-Quad: 20 faces → 12 quads (pairs of triangles).

## Compatibility Testing

Output is tested against dggridR for compatibility. Test examples:

``` r

library(dggridR)
dggs <- dgconstruct(res = 5, aperture = 3)
ref <- dgGEO_to_SEQNUM(dggs, lon, lat)$seqnum  # Reference for comparison
```

## Key Constants

Default ISEA3H orientation: - Vertex 0: (11.25°, 58.28252559°) -
Azimuth: 0°

Cell count formulas: - Aperture 3: 10 × 3^res + 2 - Aperture 4: 10 ×
4^res + 2 - Aperture 7: 10 × 7^res + 2
