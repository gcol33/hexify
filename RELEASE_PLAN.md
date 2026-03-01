# Release Plan: v0.7.0 → v1.0.0

All work happens on `experimental` branch. Each release is tested there, then merged to `main`.

---

## v0.7.0 — Spatial Analysis (bug fix + core features)

### Bug fix (critical)
- [x] Fix ap7 encoding: Class III substrate quantization instead of direct surrogate snap (was 3.6× too many cells, wrong neighbors)
- [x] Update `max_cell_id()` for ap7 bounding box
- [x] Regenerate test cache

### Features
- [x] `get_neighbors()` — k-ring neighbors for all ISEA apertures + H3
- [x] `hex_summarize()` — cell-level aggregation with tidyeval

### Status: READY — all tests passing (5417 pass, 0 fail)

---

## v0.7.1 — Patch (reserved)

Reserve for any neighbor/encoding regressions found during real-world testing of v0.7.0.

---

## v0.8.0 — Raster Integration

### Features
- [ ] `hex_extract()` — point sampling at cell centers (terra)
- [ ] `hex_zonal()` — polygon-based zonal statistics

### Before release
- [ ] Add parametrized tests for all `hex_zonal()` aggregation functions (sum, min, max, sd, count)
- [ ] Test multiple raster layers
- [ ] Test boundary polygon clipping
- [ ] Test explicit cell_id input

---

## v0.9.0 — Multi-Resolution & Topology

### Features
- [ ] `hex_compact()` / `hex_uncompact()` — lossless multi-res compression (H3 + ISEA ap7)
- [ ] `is_pentagon()` — detect 12 pentagon cells
- [ ] `hex_distance()` — grid distance between cells

### Before release
- [ ] Clarify mixed ap4/3 scope in NEWS (C++ exists but `hex_compact()` rejects non-ap7 — either expose or remove claim)
- [ ] Document `hex_distance()` NA return for cross-quad cells exceeding search depth
- [ ] Consider caching pentagon IDs per grid instead of recomputing via lonlat_to_cell loop

---

## v1.0.0 — Polish & Theory

### Features
- [ ] `hex_browse()` — interactive leaflet maps
- [ ] Theory vignette (comprehensive, 800+ lines)
- [ ] Rebuild pkgdown site

### Before release
- [ ] Fix function references in theory vignette (hexify_index_to_cell signature, etc.)
- [ ] Final R CMD check with --as-cran
- [ ] Full round of manual testing across all apertures + H3

---

## Summary

| Release | Priority | Content | Blocker |
|---|---|---|---|
| **0.7.0** | **BUG FIX** | ap7 encoding fix + neighbors + summarize | None — ready |
| 0.7.1 | Patch | Reserved for regressions | — |
| 0.8.0 | Feature | Raster extract + zonal | Need more tests |
| 0.9.0 | Feature | Compact + topology + distance | Clarify ap4/3 |
| 1.0.0 | Polish | Interactive maps + theory vignette + site | Vignette review |
