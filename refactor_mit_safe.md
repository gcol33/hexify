# Refactor Plan: MIT License Safety

This plan addresses code that references or derives from DGGRID (AGPL-3.0) to ensure hexify can be safely released under MIT license.

## Guiding Principles

1. **Mathematical formulas are not copyrightable** - but their specific implementation can be
2. **Comments are evidence** - "Based on DGGRID" comments are self-incriminating
3. **Clean-room approach** - derive from published papers, not from reading DGGRID source
4. **Document derivations** - show work from first principles to demonstrate independence

## Reference Papers (for clean-room derivation)

- Sahr, K., White, D., & Kimerling, A. J. (2003). "Geodesic Discrete Global Grid Systems"
- Snyder, J. P. (1992). "An Equal-Area Map Projection For Polyhedral Globes"
- H3 documentation (Uber, Apache 2.0 license) - similar IJK coordinate system
- Carr, D. B., Olsen, A. R., & White, D. (1992). "Hexagon Mosaic Maps for Display of Univariate and Bivariate Geographical Data"

---

## Phase 1: High Priority (Likely Derivative Code)

### Task 1.1: Rewrite index_z7.cpp and index_z7.h
**Risk: HIGH** - Header explicitly says "Exact replication of DGGRID"

**Files:**
- `src/index_z7.h` - line 3: "Exact replication of DGGRID's DgZ7StringRF.cpp logic"
- `src/index_z7.cpp` - line 3: "EXACT port of DGGRID's DgZ7StringRF.cpp"

**Action:**
1. Remove all DGGRID references from comments
2. Rewrite header comment to describe the algorithm mathematically
3. The algorithm itself (aperture-7 hierarchical indexing) is mathematical - document the derivation:
   - Aperture-7 subdivides hex into 7 children
   - Children positioned at center + 6 surrounding
   - Index encodes path through hierarchy
4. The adjacency tables (`adjacentBaseCellTable`, `inverseAdjacentBaseCellTable`) encode icosahedral topology - derive from vertex/face definitions
5. Add citation to Sahr et al. (2003) for the ISEA framework

**New header comment:**
```cpp
// index_z7.h
// Z7 Hierarchical Index for Aperture-7 Hexagonal Grids
//
// Implements hierarchical space-filling indexing for aperture-7 subdivision.
// Each parent hex divides into 7 children (1 center + 6 surrounding).
// The index encodes the traversal path from root to cell.
//
// Mathematical basis: Sahr et al. (2003) "Geodesic Discrete Global Grid Systems"
// Coordinate system: Cube coordinates (i,j,k) with constraint i+j+k=0
//
// Copyright (c) 2024-2025 hexify authors. MIT License.
```

---

### Task 1.2: Rewrite ijk_coordinates.cpp and ijk_coordinates.h
**Risk: HIGH** - "Based on DGGRID's DgIVec3D implementation"

**Files:**
- `src/ijk_coordinates.h` - line 3
- `src/ijk_coordinates.cpp` - line 3

**Action:**
1. The IJK/cube coordinate system is well-documented in multiple sources:
   - Red Blob Games: https://www.redblobgames.com/grids/hexagons/
   - H3 (Uber): https://h3geo.org/docs/core-library/coordsystems
2. Rewrite comments to cite these public sources instead of DGGRID
3. The matrix operations (`upAp7`, `downAp7`, etc.) are standard aperture-7 scaling:
   - `upAp7`: Scale down by factor of 7 (coarsen)
   - `downAp7`: Scale up by factor of 7 (refine)
   - These use the standard rotation matrices for aperture-7

**Mathematical derivation for upAp7:**
```
Aperture-7 scaling rotates by arctan(sqrt(3)/5) ≈ 19.1° and scales by sqrt(7)
Matrix: [3, -1; 1, 2] / 7 (for (i,j) -> coarser (i',j'))
```

**New header comment:**
```cpp
// ijk_coordinates.h
// Cube Coordinates for Hexagonal Grids
//
// Implements the (i,j,k) cube coordinate system for hexagonal grids
// where i + j + k = 0. This representation simplifies many hex operations.
//
// References:
// - Red Blob Games "Hexagonal Grids" (cube coordinates section)
// - H3 Coordinate Systems documentation
// - Sahr et al. (2003) for ISEA application
//
// Copyright (c) 2024-2025 hexify authors. MIT License.
```

---

### Task 1.3: Rewrite index_zorder.cpp and index_zorder.h
**Risk: MEDIUM** - "Based on DGGRID's DgZOrderStringRF.cpp"

**Files:**
- `src/index_zorder.h` - line 3
- `src/index_zorder.cpp` - line 3

**Action:**
1. Z-order (Morton) curves are a standard space-filling curve technique (1966, G.M. Morton)
2. The implementation is straightforward bit interleaving
3. Remove DGGRID reference, cite Morton's original work or standard CS textbooks

**New header comment:**
```cpp
// index_zorder.h
// Z-Order (Morton) Curve Indexing for Hexagonal Grids
//
// Implements Morton's space-filling curve (1966) for multi-resolution indexing.
// Interleaves coordinate digits to create a 1D index preserving 2D locality.
//
// Reference: Morton, G.M. (1966) "A computer oriented geodetic data base"
//
// Copyright (c) 2024-2025 hexify authors. MIT License.
```

---

### Task 1.4: Rewrite index_z3.cpp and index_z3.h
**Risk: MEDIUM** - References DGGRID function names

**Files:**
- `src/index_z3.h` - lines 3, 13, 17
- `src/index_z3.cpp` - already has good mathematical documentation

**Action:**
1. The .cpp file already has good mathematical derivation (lines 28-104)
2. Remove DGGRID function name references from header
3. The lookup table values are mathematical - derived from aperture-3 geometry
4. Keep the mathematical comments, remove DGGRID references

**New header comment:**
```cpp
// index_z3.h
// Z3 Space-Filling Index for Aperture-3 Hexagonal Grids
//
// Encodes (i,j) coordinates as a hierarchical index for aperture-3 subdivision.
// The bijective mapping preserves locality in the hex grid.
//
// Mathematical basis: Each aperture-3 level subdivides a hex into 3 children
// arranged in a triangular pattern. The encoding maps base-3 digit pairs
// to child positions using the geometric relationship.
//
// Copyright (c) 2024-2025 hexify authors. MIT License.
```

---

## Phase 2: Medium Priority (Formula References)

### Task 2.1: Clean up rcpp_seqnum.cpp comments
**Risk: MEDIUM** - Multiple "Uses DGGRID's formula" comments

**Lines to modify:**
- Line 305-309: "pattern from DGGRID"
- Line 341: "used by DGGRID's seqnum formula"
- Line 392: "Uses DGGRID's DgBoundedHexC2RF2D formula"
- Line 416: "Based on DGGRID's DgBoundedHexC3RF2D::seqNumAddress"
- Line 433: "Based on DGGRID's DgBoundedHexC3RF2D::addFromSeqNum"
- Line 561, 1471: "Calculate seqnum using DGGRID formula"
- Line 768: "Equivalent to dggridR's dgSEQNUM_to_Q2DD()"
- Line 849: "Equivalent to dggridR's dgQ2DD_to_SEQNUM()"
- Line 1077: "match dggridR vertex winding"

**Action:**
1. The seqnum formulas (`i * dim / 3`, switch on `i % 3`) are mathematical necessity:
   - For offset grids, only 1/3 (or 1/7) of coordinates are valid
   - The modulo arithmetic encodes which cells are valid
   - This is derivable from hex grid geometry
2. Replace "DGGRID's formula" with mathematical derivation explanation
3. Replace "Equivalent to dggridR" with "Produces compatible cell IDs"

**Example replacement for line 392:**
```cpp
// 2D seqnum for offset grid (aperture 3 odd resolutions)
// Only 1/3 of cells valid - those where (i+j) % 3 == 0
//
// Derivation: In Class II grids, valid cells form a triangular sublattice.
// The formula compacts the sparse 2D grid into dense 1D numbering.
// For row i, the valid j values follow pattern: j ≡ {0, 2, 1}[i % 3] (mod 3)
```

---

### Task 2.2: Clean up coordinate_transforms.cpp comments
**Risk: LOW** - vertTable is documented as "Derived from First Principles"

**Lines:**
- Line 612: "DGGRID-compatible vertTable - Derived from First Principles"
- Line 986: "Based on DGGRID's DgQuadXYtoVertex2DDConverter::compute_subtriangle"
- Line 1039, 1047: DGGRID transformation references

**Action:**
1. The vertTable derivation (lines 612-772) is excellent - keep it
2. Change "DGGRID-compatible" to "Standard ISEA"
3. Line 986: Replace with mathematical description of subtriangle classification
4. Remove specific DGGRID class/function name references

---

## Phase 3: Low Priority (Comment Cleanup)

### Task 3.1: Clean up aperture_sequence.cpp
**File:** `src/aperture_sequence.cpp`
- Line 4: "Used by DGGRID's ISEA43H grid type"

**Action:** Replace with "Implements mixed aperture grids (e.g., ISEA43H)"

---

### Task 3.2: Clean up cell_index.cpp
**File:** `src/cell_index.cpp`
- Line 2: "Unified indexing for DGGRID aperture 3, 4, 7"

**Action:** Replace with "Unified indexing for ISEA apertures 3, 4, 7"

---

### Task 3.3: Clean up sequential_numbering.h
**File:** `src/sequential_numbering.h`
- Line 9: "This provides PLANE addressing mode compatible with DGGRID"

**Action:** Replace with "Provides sequential cell numbering for ISEA grids"

---

## Phase 4: Verification

### Task 4.1: Run tests to ensure nothing breaks
```bash
Rscript -e "devtools::test()"
```

### Task 4.2: Grep for remaining DGGRID/dggridR references
```bash
grep -ri "dggrid" src/ --include="*.cpp" --include="*.h"
```

Remaining references should only be:
- Function names like `_dggrid` that are part of the API (acceptable)
- Documentation explaining compatibility (acceptable)

### Task 4.3: Update CLAUDE.md
Remove or update any DGGRID validation references that suggest code derivation.

---

## Summary Checklist

| Task | File(s) | Risk | Status |
|------|---------|------|--------|
| 1.1 | index_z7.cpp/h | HIGH | ✅ |
| 1.2 | ijk_coordinates.cpp/h | HIGH | ✅ |
| 1.3 | index_zorder.cpp/h | MEDIUM | ✅ |
| 1.4 | index_z3.cpp/h | MEDIUM | ✅ |
| 2.1 | rcpp_cell.cpp (was rcpp_seqnum.cpp) | MEDIUM | ✅ |
| 2.2 | coordinate_transforms.cpp | LOW | ✅ |
| 3.1 | aperture_sequence.cpp | LOW | ✅ |
| 3.2 | cell_index.cpp | LOW | ✅ |
| 3.3 | cell_numbering.h (was sequential_numbering.h) | LOW | ✅ |
| 3.4 | constants.h | LOW | ✅ |
| 4.1 | Tests | - | ✅ (4359 pass) |
| 4.2 | Grep verification | - | ✅ |
| 4.3 | CLAUDE.md | - | ✅ |

---

## Notes

- The projection files (projection_forward.cpp, projection_inverse.cpp) are clean - they implement Snyder's published equations
- The icosahedron.cpp file is clean - standard spherical geometry
- The aperture.cpp file is clean - standard hex quantization
- Test files that compare against dggridR are fine - testing compatibility is not code derivation
