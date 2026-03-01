# Theory Notes

Collect notes for the theory pkgdown article (`vignettes/theory.Rmd`).

---

## AP7 Encoding Bug & Fix (2026-03-01)

### The two-grid architecture

Aperture 7 grids use two coordinate systems internally:

| System | Scale (relative to quad XY) | Rotation | Role |
|---|---|---|---|
| **Substrate** | sqrt(7)^res × sqrt(7) | 0° | Class III quantized — only valid ap7 positions |
| **Surrogate** | sqrt(7)^res | −19.1° | Class I hex lattice — ±1 offsets = neighbors |

The surrogate is rotated by arctan(√(3/7)) ≈ 19.1° from the substrate and scaled down by √7. Each surrogate lattice point corresponds to exactly one valid ap7 cell. The standard 6 hex offsets {(+1,0), (+1,+1), (0,+1), (−1,0), (−1,−1), (0,−1)} in surrogate space find neighbors.

### The encoding pipeline

Correct:

```
lon/lat → Snyder forward → icosa triangle → Class III quantize → SUBSTRATE (i,j)
       → substrate_to_surrogate_ap7() → SURROGATE (i,j) → cell_id
```

The Class III quantization is the critical step. It snaps the continuous point to the nearest valid ap7 cell by:

1. Rotating to surrogate frame (−19.1°)
2. Quantizing to nearest Class I hex lattice point
3. Reading back the center, rotating to substrate frame (+19.1°)
4. Scaling up by √7 and re-quantizing to substrate

This round-trip through the surrogate ensures the final substrate position corresponds to a real ap7 cell — not an arbitrary lattice point.

### The bug

The encoding skipped Class III quantization and went directly:

```
lon/lat → scale by sqrt(7)^res → rotate −19.1° → snap to nearest hex lattice point → cell_id
```

This "snap to nearest hex lattice point" assigns cell IDs to **every** lattice position in the bounding box, not just valid ap7 cells. At resolution 3:

- Expected cells per quad: **343** (= 7³)
- Actual cells per quad: **1225** (= 35², the full bounding box)
- Phantom cells: **882** per quad — positions with cell IDs but no corresponding hexagon

The ±1 surrogate offsets then hopped between adjacent lattice points (most of which were phantoms), producing "neighbors" at ~380–450 km instead of the true nearest cells at ~80–100 km.

### The fix

Two lines in `cpp_lonlat_to_cell` and `cpp_quad_xy_to_cell`:

```cpp
// Before (wrong): direct surrogate quantization — all lattice points get IDs
hexify::quad_xy_to_surrogate_ij_ap7(qx, qy, resolution, si, sj);

// After (correct): Class III quantize → substrate → surrogate
hexify::icosa_tri_to_quad_ij(..., 7, resolution, quad, sub_i, sub_j);
hexify::substrate_to_surrogate_ap7(sub_i, sub_j, resolution, si, sj);
```

### Residual boundary excess

After the fix, cells per quad at res 3 is ~376 (not exactly 343). The ~10% excess comes from substrate positions near quad edges that map to surrogate positions slightly outside the "natural" boundary. This converges: 1.10× at res 3, 1.04× at res 4, 1.01× at res 5. These boundary cells are harmless — they have valid centers and round-trip correctly, but no geographic point maps to them from other quads.

### Key insight for the theory vignette

The Class III quantization is what makes aperture 7 work. Without it, the rotated hex lattice has ~3.6× too many positions. The quantize → inverse → rotate → re-quantize round-trip is the mathematical mechanism that selects only the 1-in-7 valid positions from the substrate lattice. This is the same mechanism described by Sahr (2008) for the ISEA aperture 7 grid — the "Class III" designation specifically refers to this rotation-based thinning.
