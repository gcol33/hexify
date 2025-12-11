# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Vision

hexify is a minimal, modern R package for generating equal-area hexagonal grids on the globe. It implements the same mathematical foundation as DGGRID (ISEA projection, hierarchical apertures) but redesigned for simplicity and modern R workflows.

**Core use case:** Given a data.frame with lon/lat columns, return the same data with:
- `hex_id` – stable DGGS cell ID
- `hex_cen_lon`, `hex_cen_lat` – cell center coordinates

No grids, no polygon generation, no complex joins. Just: "Here is my point; tell me which hex it belongs to."

## Build & Development Commands

**Important:** Use `Rscript` for command-line R execution (not `R -e`).

```bash
# Load for development
Rscript -e "devtools::load_all()"

# Run all tests
Rscript -e "devtools::test()"

# Run single test file
Rscript -e "testthat::test_file('tests/testthat/test-core_icosa.R')"

# Rebuild documentation
Rscript -e "devtools::document()"

# Compile C++ code only
Rscript -e "Rcpp::compileAttributes()"
```

## DGGRID Coordinate Systems

**Critical:** DGGRID uses multiple coordinate representations. Understanding these is essential:

### 1. PROJTRI (Triangle coordinates)
- `tnum`: Triangle number (0-19), one of 20 icosahedral faces
- `tx, ty`: Projected coordinates within that triangle (range ~[0, 1])
- This is what hexify's Snyder projection produces

### 2. Q2DI (Quad integer coordinates)
- `quad`: Quad number (0-11), pairs of triangles forming diamond shapes
- `i, j`: Integer cell indices within the quad at a given resolution
- **This is what dggridR returns for cell assignment**

### 3. SEQNUM (Sequential number)
- Single integer uniquely identifying each cell globally
- Useful for database storage

### Triangle to Quad Mapping
The 20 triangular faces are grouped into 12 quads. Each quad contains 2 triangles:
- Quad 0: North polar region
- Quads 1-10: Equatorial band
- Quad 11: South polar region

**Current issue:** hexify works in triangle space (tnum, tx, ty) but doesn't properly convert to quad space (quad, i, j) that dggridR uses.

## Current Implementation Status

### Working (verified against dggridR):
- Face detection (`which_face`) - returns correct triangle number (tnum)
- Forward Snyder projection - tx/ty match dggridR to ~1e-9 precision
- Inverse Snyder projection - lon/lat roundtrip works

### Not working:
- Cell quantization returns (face, i, j) in triangle space, not (quad, i, j) in quad space
- Cell IDs don't match dggridR's SEQNUM
- Z7 indexing has errors

## Reference Implementations

The `references/` folder contains known-good implementations:
- `dggridR-master/` - R wrapper for DGGRID
- `DGGRID-master/` - Original C++ DGGRID source

### Building DGGRID on Windows

From Git Bash:
```bash
export PATH="/c/msys64/mingw64/bin:$PATH"
cd "c:\Users\Gilles Colling\Documents\dev\hexify\references\DGGRID-master"
cmake -B build -G Ninja
cmake --build build
```
Executable: `build/src/apps/dggrid/dggrid.exe`

**Validation requirement:** ALL code must be validated against dggridR output. Do not trust existing tests - they may encode incorrect understanding. Always:
1. Generate reference data from dggridR
2. Compare hexify output to that reference
3. Only then write/update tests

Example validation:
```r
library(dggridR)
dggs <- dgconstruct(res = 5, aperture = 3)
ref <- dgGEO_to_Q2DI(dggs, lon, lat)  # Ground truth
```

## Architecture

### C++ Layer (src/)
- `core_icosa.cpp` - Icosahedron geometry, face centers, `which_face()`
- `snyder_forward.cpp` - Forward projection: lon/lat → (tnum, tx, ty)
- `snyder_inverse.cpp` - Inverse projection: (face, x, y) → lon/lat
- `hex_ap3.cpp` - Aperture-3 quantization (needs quad conversion)
- `hex_ap4.cpp` - Aperture-4 quantization
- `hex_ap7.cpp` - Aperture-7 quantization
- `hex_index*.cpp` - Space-filling curve indices

### R Layer (R/)
- `hex_assign.R` - High-level point assignment (target API)
- `dggrid.R` - DGGRID-compatible functions
- `snyder_*.R` - Projection wrappers

## Key Constants

Default ISEA3H orientation (matches dggridR):
- Vertex 0: lon=11.25°, lat=58.28252559°
- Azimuth: 0°

## Dependencies

- **LinkingTo:** Rcpp
- **Imports:** sf, Rcpp
- **Suggests:** testthat (>= 3.0.0)
