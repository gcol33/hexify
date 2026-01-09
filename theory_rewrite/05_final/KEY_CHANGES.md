# Key Changes to Aperture Diagrams Chunk

## Critical Fix: Aperture 7 Rotation

### BEFORE (Incorrect)
```r
# Center hexagon - NO ROTATION
child_r7 <- 0.52
h <- hex_vertices(0, 0, child_r7, flat_top = TRUE)  # ❌ flat_top = TRUE means 0° rotation

# Ring hexagons - NO ROTATION
for (i in 1:6) {
  ang <- pi/2 + (i-1) * pi/3
  cx <- ring_r * cos(ang)
  cy <- ring_r * sin(ang)
  h <- hex_vertices(cx, cy, child_r7, flat_top = TRUE)  # ❌ flat_top = TRUE means 0° rotation
  polygon(h$x, h$y, ...)
}
```

### AFTER (Correct - Reference-Based)
```r
# Center hexagon - ROTATED by arctan(sqrt(3/7)) ≈ 19.1°
child_r7 <- 1.7 / sqrt(7)  # ✓ Correct scale factor
rot_angle <- atan(sqrt(3/7))  # ✓ 19.1° rotation from Sahr (2008)
angles_center <- seq(rot_angle, 2*pi + rot_angle, length.out = 7)
h_center <- list(x = child_r7 * cos(angles_center),
                 y = child_r7 * sin(angles_center))
polygon(h_center$x, h_center$y, ...)

# Ring hexagons - SAME ROTATION as center
ring_r <- 2 * child_r7 * cos(pi/6)  # ✓ Correct geometric distance
for (i in 1:6) {
  ang <- pi/2 + (i-1) * pi/3
  cx <- ring_r * cos(ang)
  cy <- ring_r * sin(ang)

  # ✓ Apply same rotation to all ring hexagons
  angles_ring <- seq(rot_angle, 2*pi + rot_angle, length.out = 7)
  h_ring <- list(x = cx + child_r7 * cos(angles_ring),
                 y = cy + child_r7 * sin(angles_ring))
  polygon(h_ring$x, h_ring$y, ...)
}
```

## Scale Factor Corrections

| Aperture | Before (Hardcoded) | After (Correct Formula) | Reference |
|----------|-------------------|-------------------------|-----------|
| 3 | `child_r = 0.85` | `child_r = 1.7 / sqrt(3)` ≈ 0.98 | Linear scale √3 |
| 4 | `child_r4 = 0.7` | `child_r4 = 1.7 / 2` = 0.85 | Linear scale 2.0 |
| 7 | `child_r7 = 0.52` | `child_r7 = 1.7 / sqrt(7)` ≈ 0.64 | Linear scale √7 |

## Why This Matters

**Aperture 7 is the most significant fix:**
- The previous code showed NO rotation (flat_top = TRUE)
- The reference states rotation = arctan(√(3/7)) ≈ 19.1°
- This rotation is a **defining characteristic** of aperture 7 grids
- Without it, the diagram misrepresents the actual grid structure

**From Sahr (2008, p. 176):**
> "For aperture 7, each level adds a rotation of arctan(√(3/7)) ≈ 19.1°"

The visualization now correctly demonstrates this critical geometric property.

## Visual Impact

### Panel 1 (Aperture 3)
- Child hexagons slightly larger and better positioned
- 30° rotation correctly represented (pointy-top orientation)

### Panel 2 (Aperture 4)
- Child hexagons slightly larger and better aligned
- No rotation correctly maintained (flat-top orientation)

### Panel 3 (Aperture 7) - **MAJOR CHANGE**
- Child hexagons larger and correctly sized
- **All 7 hexagons now show 19.1° rotation** (previously showed 0°)
- Ring distance geometrically derived instead of hardcoded
- Accurately represents the "Class III" skewed orientation
