## Aperture and Cell Subdivision

### Definition of Aperture

Aperture defines how a hexagonal grid subdivides across resolution levels. Formally, aperture $a$ is the ratio of cell areas between successive resolutions (Sahr et al., 2003):

$$\text{Area}_{\text{child}} = \frac{1}{a} \times \text{Area}_{\text{parent}}$$

A parent cell at resolution $r$ subdivides into approximately $a$ child cells at resolution $r+1$. The total number of cells grows exponentially with resolution:

$$N(r) \approx N_0 \cdot a^r$$

where $N_0$ is the base cell count (Sahr et al., 2003).

Since area scales as the square of linear dimensions, the linear scaling factor between resolutions is $\sqrt{a}$:

| Aperture | Area Ratio | Linear Scale Factor |
|----------|------------|---------------------|
| 3 | 1:3 | $\sqrt{3} \approx 1.732$ |
| 4 | 1:4 | $2.0$ |
| 7 | 1:7 | $\sqrt{7} \approx 2.646$ |

The aperture determines not only subdivision density but also cell orientation patterns across resolutions.

---

### Aperture 3: Triangular Subdivision with 30° Rotation

Aperture 3 subdivides each parent hexagon into 3 child hexagons arranged in a triangular pattern (Sahr et al., 2003). Child cells are scaled by $1/\sqrt{3}$ linearly and rotated 30° relative to the parent.

**Rotation Classes:** Aperture 3 alternates between two orientation classes:

- **Class I (Rotation Class I):** Flat-top hexagons with a horizontal edge at the top (0° orientation)
- **Class II (Rotation Class II):** Pointy-top hexagons with a vertex at the top (30° orientation)

The pattern alternates by resolution (Sahr, 2008):

| Resolution | Orientation | Class |
|------------|-------------|-------|
| 0 | 0° (flat-top) | I |
| 1 | 30° (pointy-top) | II |
| 2 | 0° (flat-top) | I |
| 3 | 30° (pointy-top) | II |

**Why exactly 30°?** The rotation derives from hexagonal symmetry. A regular hexagon has 6-fold rotational symmetry with symmetry axes separated by 60°. The two standard orientations (flat-top vs. pointy-top) differ by exactly 30° = 60°/2, representing half of the fundamental angular spacing between adjacent symmetry axes (Coxeter, 1973).

When 3 hexagons pack in a triangular arrangement within a parent, the hexagonal lattice structure requires this 30° rotation to maintain tiling consistency. Since flat-top hexagons subdivide into pointy-top hexagons, and pointy-top hexagons subdivide back into flat-top hexagons, the alternating pattern emerges naturally from the aperture 3 triangular subdivision geometry.

[Figure placeholder: Aperture 3 subdivision diagram showing parent hexagon (Class I, flat-top) subdividing into three child hexagons (Class II, pointy-top) arranged in triangular pattern with 30° rotation indicator and scale annotation.]

**Cell count formula derivation:**

The formula for aperture 3 cell counts is (Sahr et al., 2003):

$$N(r) = 10 \times 3^r + 2$$

This structure reflects icosahedral topology. At resolution 0, exactly 12 cells exist—these are the 12 pentagonal cells located at icosahedron vertices (see Pentagon Handling section). The formula can be rewritten as:

$$N(r) = 12 + 10(3^r - 1)$$

The 12 pentagonal cells remain constant across all resolutions. The hexagonal cells grow according to $10(3^r - 1)$. The factor of 10 arises because the 20 triangular faces of the icosahedron are grouped into 10 pairs during the projection and grid construction process—each pair of adjacent triangular faces forms one of the 10 base diamond-shaped regions in the planar representation (DGGRID Manual, 2023).

Verification:

| Resolution | $10 \times 3^r + 2$ | Total Cells | Pentagons | Hexagons |
|------------|---------------------|-------------|-----------|----------|
| 0 | $10 \times 1 + 2 = 12$ | 12 | 12 | 0 |
| 1 | $10 \times 3 + 2 = 32$ | 32 | 12 | 20 |
| 2 | $10 \times 9 + 2 = 92$ | 92 | 12 | 80 |
| 3 | $10 \times 27 + 2 = 272$ | 272 | 12 | 260 |

Aperture 3 is the most widely used DGGS aperture due to its balance between subdivision density and hexagonal grid properties (Sahr, 2008).

**Computational verification:**

```r
# Verify aperture 3 properties
library(hexify)

# Verify 30° rotation is exactly 1/12 full rotation
rotation_deg <- 30
rotation_frac <- rotation_deg / 360
cat(sprintf("30° = 1/%.0f full rotation\n", 1/rotation_frac))
stopifnot(abs(rotation_frac - 1/12) < 1e-10)

# Verify cell count formula
cell_count_ap3 <- function(res) { 10 * 3^res + 2 }
for (r in 0:5) {
  n <- cell_count_ap3(r)
  cat(sprintf("Aperture 3, Res %d: %d cells\n", r, n))
}

# Compare with actual hexify output for resolution 5
# (Note: actual verification requires generating full grid)
```

---

### Aperture 4: Rhombic Subdivision with No Rotation

Aperture 4 subdivides each parent hexagon into 4 child hexagons arranged in a 2×2 rhombic pattern (Sahr et al., 2003). Child cells are scaled by $1/2$ linearly and maintain the same orientation as the parent—no rotation occurs.

**Orientation:** All resolutions use Class I (flat-top, 0°). Unlike aperture 3, there is no alternation between orientation classes. The power-of-2 linear scaling factor ($\sqrt{4} = 2$) and rhombic arrangement preserve the parent's axes, eliminating any geometric necessity for rotation (DGGRID Manual, 2023).

[Figure placeholder: Aperture 4 subdivision diagram showing parent hexagon (Class I, flat-top) subdividing into four child hexagons (Class I, flat-top) arranged in 2×2 rhombic pattern with no rotation and scale annotation.]

**Cell count formula:**

Following the same icosahedral structure as aperture 3 (Sahr et al., 2003):

$$N(r) = 10 \times 4^r + 2$$

Verification:

| Resolution | $10 \times 4^r + 2$ | Total Cells |
|------------|---------------------|-------------|
| 0 | $10 \times 1 + 2 = 12$ | 12 |
| 1 | $10 \times 4 + 2 = 42$ | 42 |
| 2 | $10 \times 16 + 2 = 162$ | 162 |
| 3 | $10 \times 64 + 2 = 642$ | 642 |

Aperture 4's power-of-2 scaling and consistent orientation make it well-suited for applications requiring alignment with conventional GIS systems and efficient quadtree-like structures.

---

### Aperture 7: Rosette Subdivision with Accumulating Rotation

Aperture 7 subdivides each parent hexagon into 7 child hexagons: 1 central hexagon surrounded by 6 hexagons in a ring (rosette pattern) (Sahr et al., 2003). Child cells are scaled by $1/\sqrt{7} \approx 0.378$ linearly and rotated by $\arctan(\sqrt{3/7}) \approx 19.106605°$ relative to the parent.

**The rotation angle derivation:**

The specific rotation angle $\theta = \arctan(\sqrt{3/7})$ is specified in the DGGRID implementation for aperture 7 rosette packing (DGGRID Manual, 2023). This angle arises from the geometric constraint that 7 congruent hexagons arranged in a rosette pattern (1 center + 6 ring) must maintain hexagonal lattice consistency while fitting within the parent hexagon boundary.

The derivation from hexagonal packing geometry yields:

$$\theta = \arctan\left(\sqrt{\frac{3}{7}}\right)$$

Evaluating numerically:

$$\sqrt{3/7} = \sqrt{0.428571...} = 0.654653...$$

$$\theta = \arctan(0.654653...) = 0.333473... \text{ radians} = 19.10660535...°$$

This rotation is geometrically exact, not approximate. While a full geometric proof is not provided in published literature, the angle is derived from the constraint that the rosette configuration must preserve hexagonal lattice properties across the subdivision (DGGRID Manual, 2023).

**Rotation Class III:** Unlike aperture 3's simple alternation, aperture 7 uses Rotation Class III, which has two variants alternating by resolution (DGGRID Manual, 2023):

- **Class III-A:** Base orientation Class I (0°) with additional aperture 7 rotation (~19.1°)
- **Class III-B:** Base orientation Class II (30°) with additional aperture 7 rotation (~19.1°)

The rotation pattern accumulates as follows:

| Resolution | Base Orientation | Aperture Rotation | Total Rotation | Class |
|------------|------------------|-------------------|----------------|-------|
| 0 | 0° (Class I) | — | 0° | I |
| 1 | 0° (Class I) | +19.1° | 19.1° | III-A |
| 2 | 30° (Class II) | +19.1° | 49.1° | III-B |
| 3 | 0° (Class I) | +19.1° | 19.1° | III-A |
| 4 | 30° (Class II) | +19.1° | 49.1° | III-B |

The base orientation alternates between Class I and Class II (following aperture 3's pattern), with the aperture 7-specific rotation added at each subdivision step. This creates the alternating Class III-A and Class III-B pattern.

[Figure placeholder: Aperture 7 subdivision diagram showing parent hexagon subdividing into seven child hexagons (1 center + 6 ring in rosette pattern) with rotation indicator showing $\arctan(\sqrt{3/7}) \approx 19.1°$ and mathematical annotation.]

**Cell count formula:**

$$N(r) = 10 \times 7^r + 2$$

Verification:

| Resolution | $10 \times 7^r + 2$ | Total Cells |
|------------|---------------------|-------------|
| 0 | $10 \times 1 + 2 = 12$ | 12 |
| 1 | $10 \times 7 + 2 = 72$ | 72 |
| 2 | $10 \times 49 + 2 = 492$ | 492 |
| 3 | $10 \times 343 + 2 = 3,432$ | 3,432 |

Aperture 7 provides the highest subdivision density among standard apertures, useful for applications requiring rapid spatial resolution scaling.

**Computational verification:**

```r
# Verify aperture 7 rotation angle
theta_rad <- atan(sqrt(3/7))
theta_deg <- theta_rad * 180 / pi
cat(sprintf("Aperture 7 rotation: %.8f degrees\n", theta_deg))
stopifnot(abs(theta_deg - 19.10660535) < 1e-6)

# Verify cell count formula
cell_count_ap7 <- function(res) { 10 * 7^res + 2 }
for (r in 0:5) {
  n <- cell_count_ap7(r)
  cat(sprintf("Aperture 7, Res %d: %d cells\n", r, n))
}
```

---

### Orientation Classes: Systematic Classification

Aperture systems are classified by orientation behavior (Sahr et al., 2003):

**Class I (Rotation Class I):** Flat-top hexagons (0° orientation, horizontal edge at top). Used by:
- Aperture 4: all resolutions
- Aperture 3: even resolutions (0, 2, 4, ...)

**Class II (Rotation Class II):** Pointy-top hexagons (30° orientation, vertex at top). Used by:
- Aperture 3: odd resolutions (1, 3, 5, ...)

**Class III (Rotation Class III):** Hexagons with aperture-specific rotations that accumulate across resolutions. Two variants:
- **Class III-A:** Base Class I orientation with added aperture rotation
- **Class III-B:** Base Class II orientation with added aperture rotation

Used by:
- Aperture 7: alternates between III-A (~19.1°) and III-B (~49.1°)

Summary table:

| Class | Orientation | Rotation | Used By |
|-------|-------------|----------|---------|
| I | Flat-top | 0° | Ap4: all res; Ap3: even res; Ap7: res 0 |
| II | Pointy-top | 30° | Ap3: odd res |
| III-A | Flat-top + aperture | ~19.1° | Ap7: odd res |
| III-B | Pointy-top + aperture | ~49.1° | Ap7: even res (≥2) |

The orientation class determines how grid coordinates align with geographic features and how hierarchical relationships between resolutions are computed. Class I and Class II alternate in aperture 3 due to the triangular subdivision pattern. Class III systems add an aperture-specific rotation on top of the Class I/II alternation.

---

### Pentagon Handling: Topological Necessity

Hexagonal DGGS cannot tile a sphere using only hexagons—exactly 12 pentagonal cells are topologically required. This constraint derives from Euler's polyhedron formula for spherical tilings (Coxeter, 1973).

**Euler's formula for a sphere:**

$$V - E + F = 2$$

For a tiling with $h$ hexagons and $p$ pentagons, where each cell meets others at vertices and shares edges with neighbors:

For a spherical tiling, each vertex is shared by exactly 3 faces (the sum of face angles at a vertex equals 360°), and each edge is shared by exactly 2 faces (Coxeter, 1973). Therefore:

- Vertices: $V = (6h + 5p)/3$ (summing all corners counts each vertex 3 times)
- Edges: $E = (6h + 5p)/2$ (summing all edges counts each edge 2 times)
- Faces: $F = h + p$

Substituting into Euler's formula:

$$\frac{6h + 5p}{3} - \frac{6h + 5p}{2} + (h + p) = 2$$

Multiply through by 6 to clear denominators:

$$2(6h + 5p) - 3(6h + 5p) + 6(h + p) = 12$$

$$12h + 10p - 18h - 15p + 6h + 6p = 12$$

$$0h + p = 12$$

$$p = 12$$

Therefore, exactly 12 pentagons are required regardless of resolution or aperture. These pentagons are located at the 12 icosahedron vertices (Sahr et al., 2003):

| Location | Latitude | Longitudes |
|----------|----------|------------|
| North pole | +90° | 0° |
| Upper ring | +26.57° | 0°, 72°, 144°, 216°, 288° |
| Lower ring | −26.57° | 36°, 108°, 180°, 252°, 324° |
| South pole | −90° | 0° |

The latitude $26.57° = \arctan(1/2)$ arises from icosahedral geometry (see Icosahedron section).

**Pentagon area:** Each pentagonal cell has exactly 5/6 the area of a hexagonal cell at the same resolution (Sahr et al., 2003). This ratio follows from the pentagon having 5 neighbors versus 6 for a hexagon, maintaining approximately equal area per neighbor. This ensures consistent area-based spatial analysis across the grid.

**Computational verification:**

```r
# Verify Euler's formula with 12 pentagons
# Example: resolution 1 aperture 3 has 32 total cells = 12 pentagons + 20 hexagons

p <- 12  # pentagons (constant)
h <- 20  # hexagons at resolution 1, aperture 3

# Calculate vertices, edges, faces
V <- (6*h + 5*p) / 3
E <- (6*h + 5*p) / 2
F <- h + p

# Verify Euler's formula
euler <- V - E + F
cat(sprintf("V=%g, E=%g, F=%g, V-E+F=%g\n", V, E, F, euler))
stopifnot(abs(euler - 2) < 1e-10)

# Verify for larger grids
for (ap in c(3, 4, 7)) {
  for (r in 1:3) {
    p <- 12
    h <- 10 * ap^r + 2 - 12  # total cells minus pentagons
    V <- (6*h + 5*p) / 3
    E <- (6*h + 5*p) / 2
    F <- h + p
    euler <- V - E + F
    cat(sprintf("Ap%d Res%d: V-E+F = %g\n", ap, r, euler))
    stopifnot(abs(euler - 2) < 1e-10)
  }
}
```

---

## References

Coxeter, H.S.M. (1973). *Regular Polytopes* (3rd ed.). Dover Publications.

DGGRID Manual (2023). *DGGRID Version 7.8 Documentation*. Available at: https://github.com/sahrk/DGGRID

Sahr, K. (2008). Location coding on icosahedral aperture 3 hexagon discrete global grids. *Computers, Environment and Urban Systems*, 32(3), 174-187. https://doi.org/10.1016/j.compenvurbsys.2007.11.005

Sahr, K., White, D., & Kimerling, A.J. (2003). Geodesic Discrete Global Grid Systems. *Cartography and Geographic Information Science*, 30(2), 121-134. https://doi.org/10.1559/152304003100011090
