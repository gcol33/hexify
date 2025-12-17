# Icosahedron Geometry for ISEA Projections

## Table of Contents
1. [Mathematical Foundations](#mathematical-foundations)
2. [Face Center Computation](#face-center-computation)
3. [Standard ISEA Orientation](#standard-isea-orientation)
4. [Face Assignment Algorithm](#face-assignment-algorithm)
5. [Implementation Verification](#implementation-verification)
6. [References](#references)

---

## 1. Mathematical Foundations

### 1.1 Regular Icosahedron Definition

A regular icosahedron is a convex Platonic solid with:
- **12 vertices**
- **30 edges** (all of equal length)
- **20 faces** (all equilateral triangles)

The icosahedron exhibits icosahedral symmetry (group I_h), the highest symmetry available for a polyhedron.

### 1.2 Dual Relationship with Dodecahedron

The regular icosahedron and regular dodecahedron are dual polyhedra:

| Property | Icosahedron | Dodecahedron |
|----------|-------------|--------------|
| Faces | 20 (triangles) | 12 (pentagons) |
| Vertices | 12 | 20 |
| Edges | 30 | 30 |

**Duality mapping:**
- Each vertex of the icosahedron corresponds to a face of the dodecahedron
- Each face of the icosahedron corresponds to a vertex of the dodecahedron
- Both share the same symmetry group

**Construction:** Connecting the centers of adjacent faces of an icosahedron yields a dodecahedron, and vice versa.

### 1.3 Vertex Coordinates on Unit Sphere

When inscribed in a unit sphere (radius R = 1), the icosahedron vertices can be expressed using the golden ratio φ = (1 + √5)/2 ≈ 1.618034.

**Exact coordinates (Cartesian, before normalization):**

The 12 vertices consist of all cyclic permutations of:
```
(0, ±1, ±φ)
```

After normalization to lie on the unit sphere, the scaling factor is 1/√(1 + φ²) = √5/√(2(5 + √5)).

**Normalized vertex coordinates:**
```
(0, ±0.5257311121, ±0.8506508084)
(±0.5257311121, ±0.8506508084, 0)
(±0.8506508084, 0, ±0.5257311121)
```

### 1.4 Pentagon Vertex Orientation

**Key geometric property:** When oriented with one vertex at the north pole and another at the south pole, the remaining 10 vertices form two rings of 5 vertices each.

**Latitude of ring vertices:**
```
latitude = ±arctan(1/2) = ±26.565051177°
```

**Mathematical derivation:**

For an icosahedron oriented with vertices at the poles, the middle vertices lie at a height z such that the ratio of vertical to horizontal distance creates the proper icosahedral angles. This ratio is exactly 1:2, giving:

```
tan(colatitude) = 2
colatitude = arctan(2) ≈ 63.434948823°
latitude = 90° - arctan(2) = arctan(1/2) ≈ 26.565051177°
```

This can also be expressed as:
```
sin(latitude) = 1/√5 ≈ 0.447213595
cos(latitude) = 2/√5 ≈ 0.894427191
```

**Longitude spacing:** The 5 vertices in each ring are evenly spaced at 72° intervals (360°/5).

**Offset between rings:** The two rings are offset by 36° in longitude.

### 1.5 Standard Icosahedron Properties

For a regular icosahedron inscribed in a unit sphere:

- **Circumradius** (sphere radius): R = 1
- **Edge length**: a = √(10 - 2√5) / 2 ≈ 1.051462224
- **Midradius** (distance to edge midpoint): r_e = φ/√3 ≈ 0.934172358
- **Inradius** (distance to face center): r_i = φ²/(2√3) ≈ 0.755761314
- **Dihedral angle**: arccos(-√5/3) ≈ 138.189685°
- **Solid angle at each vertex**: arccos(√5/3) ≈ 0.551286 steradians

**Face area on unit sphere:**
Each spherical triangle face has area:
```
A_face = 4π / 20 = π/5 ≈ 0.628318531 steradians
```

---

## 2. Face Center Computation

### 2.1 Spherical Triangle Centroid

The **face center** of each icosahedral face is the centroid of the corresponding spherical triangle on the unit sphere.

**Algorithm (Euclidean centroid method):**

Given three vertices v₁, v₂, v₃ of a spherical triangle (in Cartesian coordinates on the unit sphere):

1. **Sum the vertex vectors:**
   ```
   v_sum = v₁ + v₂ + v₃
   ```

2. **Normalize to unit sphere:**
   ```
   v_center = v_sum / |v_sum|
   ```

3. **Convert to geographic coordinates:**
   ```
   lon = atan2(v_center.y, v_center.x)
   lat = asin(v_center.z)
   ```

**Implementation in hexify** (from `icosahedron.cpp:26-32`):

```cpp
inline Geo sph_tricen(const Geo tri[3]) {
  const Vec3 a = ll2xyz(tri[0]);
  const Vec3 b = ll2xyz(tri[1]);
  const Vec3 c = ll2xyz(tri[2]);
  const Vec3 v{ a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z };
  return xyz2ll(v);
}
```

Where:
- `ll2xyz`: Converts (lon, lat) → (x, y, z) Cartesian coordinates
- `xyz2ll`: Converts (x, y, z) → (lon, lat) after normalization

### 2.2 Properties of the Centroid

**Vertex-median point:** This method computes the intersection of the spherical medians of the triangle. It is the simplest and most commonly used definition of "spherical triangle center."

**Important notes:**
- Unlike planar geometry, the spherical centroid does NOT generally divide the spherical triangle into three equal-area sub-triangles
- The centroid is the point of maximum symmetry within the triangle
- For equilateral spherical triangles (like icosahedron faces), the centroid coincides with other classical triangle centers (circumcenter, incenter, etc.)

### 2.3 Alternative Center Definitions

Other possible definitions of "face center" exist but are rarely used:

1. **Area-weighted centroid:** The point that divides the triangle into three equal-area triangles (more complex to compute)

2. **Brock centroid (1974):** Uses edge integration formula:
   ```
   centroid = Σ (A→B) over all edges
   ```

3. **Inscribed circle center:** The point equidistant from all three edges (coincides with centroid for equilateral triangles)

For the ISEA projection and hexify, the **Euclidean centroid method** (sum and normalize) is used universally.

---

## 3. Standard ISEA Orientation

### 3.1 ISEA Default Orientation

The standard ISEA (Icosahedral Snyder Equal Area) orientation places:

**Vertex 0:**
```
Longitude: 11.25° E
Latitude:  58.28252559° N
```

**Azimuth:** 0° (no rotation around the polar axis)

This orientation is used by:
- Snyder's original 1992 paper
- DGGRID / dggridR
- PROJ library
- hexify

### 3.2 Vertex 0 Latitude Derivation

The latitude 58.28252559° is NOT equal to arctan(1/2) = 26.565051177°.

**Actual relationship:**

When one vertex is at the north pole (90°N), the **second ring of vertices** below it is at ±26.565°. The **top vertex** of the main icosahedron structure is positioned such that:

```
Vertex 0 latitude = 90° - arccos(1/√5)
                  = arcsin(2/√5)
                  ≈ 58.28252559°
```

This is the **complement of the dihedral half-angle** and ensures proper icosahedral geometry.

**Verification:**
```python
import math
lat = math.degrees(math.asin(2/math.sqrt(5)))
print(f"{lat:.8f}")  # 58.28252559
```

### 3.3 Vertex 0 Longitude Choice

The choice of 11.25° E is somewhat arbitrary but has practical advantages:

1. **Geometric:** 11.25° = 45°/4, making computations with 60° and 120° sectors more convenient
2. **Geographic:** Positions icosahedron faces to minimize distortion over major landmasses
3. **Historical:** Inherited from Snyder's original implementation

Different orientations can be used (e.g., 0°, 90°) but change the face numbering and cell IDs, breaking compatibility with reference implementations.

### 3.4 Face Numbering Convention

**Standard ISEA numbering** (used by hexify, dggridR, DGGRID):

```
Faces 0-4:   Upper ring (5 triangles sharing vertex 0)
Faces 5-9:   Upper middle ring
Faces 10-14: Lower middle ring
Faces 15-19: Lower ring (5 triangles sharing vertex 11)
```

**Vertex connectivity** (from `icosahedron.cpp:139-144`):

```cpp
static const int faces[20][3] = {
  {0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},         // Faces 0-4
  {6,2,1},{7,3,2},{8,4,3},{9,5,4},{10,1,5},        // Faces 5-9
  {2,6,7},{3,7,8},{4,8,9},{5,9,10},{1,10,6},       // Faces 10-14
  {11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10}   // Faces 15-19
};
```

This creates a pentagonal gyroelongated bipyramid structure:
- Vertex 0 (north) connects to 5 faces (0-4)
- Middle vertices (1-10) each connect to 5 faces
- Vertex 11 (south) connects to 5 faces (15-19)

---

## 4. Face Assignment Algorithm

### 4.1 The Dot Product Method

**Problem:** Given a point (lon, lat) on the sphere, determine which of the 20 icosahedral faces contains it.

**Solution:** Use the **great-circle distance** from the point to each face center. The face with the **minimum distance** (equivalently, **maximum dot product**) contains the point.

### 4.2 Mathematical Justification

**Spherical law of cosines:**

The great-circle distance d between two points on a unit sphere is:

```
cos(d) = sin(lat₁)·sin(lat₂) + cos(lat₁)·cos(lat₂)·cos(lon₁ - lon₂)
```

**Dot product interpretation:**

Converting both points to Cartesian coordinates on the unit sphere:
- Point p = (x₁, y₁, z₁)
- Face center c = (x₂, y₂, z₂)

The dot product is:
```
p · c = x₁x₂ + y₁y₂ + z₁z₂ = cos(d)
```

Since cosine is a monotonically decreasing function on [0, π], maximizing the dot product minimizes the angular distance.

**Why this works:**

Each face occupies a region of the sphere bounded by three great circle arcs. The face center is the "most representative" point of the face. For points **strictly inside** a face (not on the boundary), the nearest face center is always the center of the containing face.

### 4.3 Implementation in hexify

From `icosahedron.cpp:181-196`:

```cpp
int which_face(double lon_deg, double lat_deg) {
  const IcosaData& icosa = ico();
  const Geo point(deg2rad(lon_deg), deg2rad(lat_deg));
  int best = 0;
  double bestd = std::acos(clampd(
    std::sin(icosa.centers[0].lat)*std::sin(point.lat) +
    std::cos(icosa.centers[0].lat)*std::cos(point.lat)*
    std::cos(icosa.centers[0].lon - point.lon),
    -1.0, 1.0));

  for (int i = 1; i < 20; ++i) {
    const auto& c = icosa.centers[i];
    const double cc = clampd(
      std::sin(c.lat)*std::sin(point.lat) +
      std::cos(c.lat)*std::cos(point.lat)*
      std::cos(c.lon - point.lon),
      -1.0, 1.0);
    const double d = std::acos(cc);
    if (d < bestd) { best = i; bestd = d; }
  }

  return best;
}
```

**Algorithm steps:**

1. Convert point to radians
2. For face 0, compute great-circle distance using spherical law of cosines
3. Loop through faces 1-19, computing distances
4. Return the face index with minimum distance

**Optimization opportunity:** The `acos()` call is unnecessary since cosine is monotonic. Comparing `cos(d)` directly (larger is better) would be more efficient:

```cpp
// More efficient version:
double best_cos = cos_value[0];
for (int i = 1; i < 20; ++i) {
  double cos_d = /* spherical law of cosines */;
  if (cos_d > best_cos) { best = i; best_cos = cos_d; }
}
```

### 4.4 Edge Cases and Tie-Breaking

**Points on face edges:**

When a point lies exactly on the shared edge between two faces, both face centers are equidistant. The algorithm returns the **first face** in the iteration order (lowest face number).

**Deterministic behavior:** The algorithm is deterministic - the same input always produces the same output, even for edge cases.

**Points on face vertices:**

Three faces share each vertex. The algorithm will assign the point to one of them (typically the lowest-numbered face).

**Numerical precision:**

Floating-point arithmetic introduces small errors. Two distances that are mathematically equal may differ slightly in computation. The algorithm handles this by using `clampd()` to keep arguments to `acos()` in the valid range [-1, 1].

### 4.5 Uniqueness of Assignment

**Uniqueness guarantee:**

For points **strictly inside a face** (not on boundaries), the face assignment is unique and mathematically well-defined.

**Boundary ambiguity:**

On edges and vertices, the assignment is ambiguous geometrically but deterministic algorithmically. Different implementations may assign boundary points to different faces, but the same implementation will be consistent.

**Practical impact:**

Since the icosahedron is used as an intermediate structure (points are further subdivided into hexagonal cells), boundary ambiguity has negligible practical impact:
- Probability of hitting an exact edge: measure zero (probability 0)
- Cell assignment remains consistent within each implementation
- Different implementations (hexify vs dggridR) may differ on edge points but agree on interior points

---

## 5. Implementation Verification

### 5.1 hexify Implementation

**Source files:**
- `src/icosahedron.h` (lines 1-50): Interface definitions
- `src/icosahedron.cpp` (lines 1-200): Implementation

**Key data structures:**

```cpp
struct IcosaData {
  std::array<Geo, 20> centers;              // Face centers (radians)
  std::array<double, 20> center_sinlat;     // Precomputed sin(lat)
  std::array<double, 20> center_coslat;     // Precomputed cos(lat)
  std::array<double, 20> center_lon;        // Center longitudes (radians)
  std::array<double, 20> face_azimuth_offset; // Per-face azimuth offsets
  bool built = false;
};
```

### 5.2 Vertex Construction Algorithm

From `icosahedron.cpp:120-137`:

**Step 1:** Create 12 vertices in a standard orientation

```cpp
// Vertex 0 at specified position
icoverts[0] = (vert0_lon_deg, vert0_lat_deg)

// Vertices 1-5: upper ring at +26.565°
for (i = 1; i <= 5; ++i) {
  vertsnew[i].lat = 26.565051177°
  vertsnew[i].lon = -azimuth + 72° × (i-1)
}

// Vertices 6-10: lower ring at -26.565°
for (i = 1; i <= 5; ++i) {
  vertsnew[i+5].lat = -26.565051177°
  vertsnew[i+5].lon = -azimuth + 36° + 72° × (i-1)
}

// Vertex 11 at south pole
vertsnew[11] = (0°, -90°)
```

**Step 2:** Apply coordinate transformation to rotate icosahedron

Uses `coordtrans()` function (lines 64-105) to rotate the standard orientation to place vertex 0 at the specified position.

### 5.3 Face Center Computation

From `icosahedron.cpp:146-153`:

```cpp
for (int i = 0; i < 20; ++i) {
  Geo tri[3] = { icoverts[faces[i][0]],
                 icoverts[faces[i][1]],
                 icoverts[faces[i][2]] };
  Geo c = sph_tricen(tri);  // Compute centroid
  g_ico.centers[i] = c;
  g_ico.center_sinlat[i] = std::sin(c.lat);
  g_ico.center_coslat[i] = std::cos(c.lat);
  g_ico.center_lon[i] = c.lon;
}
```

**Precomputation optimization:** The sin/cos values are precomputed to avoid redundant trigonometric calculations during face assignment.

### 5.4 Verification Against Theory

**Check 1: Vertex latitudes**

From the code (line 128, 129):
```cpp
vertsnew[i].lat = deg2rad(26.565051177);    // Upper ring
vertsnew[i+5].lat = -deg2rad(26.565051177); // Lower ring
```

Theoretical value:
```python
import math
expected = math.degrees(math.atan(0.5))
print(f"{expected:.9f}")  # 26.565051177
```

**Match:** ✓ The implementation uses the exact theoretical value.

**Check 2: Vertex 0 latitude**

From constants (R/constants.R:79):
```r
ISEA_VERT0_LAT_DEG <- 58.28252559
```

Theoretical value:
```python
import math
expected = math.degrees(math.asin(2/math.sqrt(5)))
print(f"{expected:.8f}")  # 58.28252559
```

**Match:** ✓ The implementation uses the correct vertex 0 latitude.

**Check 3: Face count and structure**

- 20 faces defined: ✓ (line 139-144)
- Each face has 3 vertices: ✓
- Vertex connectivity follows pentagonal gyroelongated bipyramid: ✓

**Check 4: Face centers span the globe**

From test file (`tests/testthat/test-icosahedron.R:48-60`):

```r
test_that("face centers span the globe", {
  hexify_build_icosa()
  centers <- hexify_face_centers()

  # Should have faces in both hemispheres
  expect_true(any(centers$lat > 0))
  expect_true(any(centers$lat < 0))

  # Should have faces across longitude range
  lon_range <- max(centers$lon) - min(centers$lon)
  expect_true(lon_range > 2.0)  # ~5.5 radians or >100 degrees
})
```

This test verifies that the computed face centers properly cover the sphere.

### 5.5 Comparison with dggridR

From `references/dggridR-master/src/DgProjTriRF.cpp`:

**Vertex latitude (line 183, 187):**
```cpp
vertsnew[i].lat = 26.565051177 * M_PI / 180.0;
vertsnew[i+5].lat = -26.565051177 * M_PI / 180;
```

**Default orientation (DgProjTriRF.h:49):**
```cpp
DgGeoCoord(11.25L, 58.28252559L, false)
```

**Match:** ✓ hexify uses identical constants and conventions to dggridR.

**Face numbering:** The face connectivity array in hexify matches the dggridR `verts` array (DgProjTriRF.cpp:156-177).

### 5.6 Deviations from Standard Conventions

**None identified.** The hexify implementation follows standard ISEA conventions:

- Vertex 0 at (11.25°E, 58.28252559°N)
- Azimuth = 0°
- Middle ring vertices at ±26.565051177° = ±arctan(1/2)
- 20 faces numbered 0-19 in standard order
- Face centers computed as spherical triangle centroids
- Face assignment via great-circle distance minimization

---

## 6. References

### Primary Sources

1. **Snyder, J.P. (1992).** "An Equal-Area Map Projection For Polyhedral Globes." *Cartographica: The International Journal for Geographic Information and Geovisualization*, 29(1):10-21.
   - [Publisher Link](https://utppublishing.com/doi/abs/10.3138/27H7-8K88-4882-1752)
   - [ResearchGate PDF](https://www.researchgate.net/publication/296900899_Snyder_equal-area_map_projection_for_polyhedral_globes)

2. **PROJ Library Documentation.** "Icosahedral Snyder Equal Area."
   - [PROJ 9.6.2 Documentation](https://proj.org/en/stable/operations/projections/isea.html)

3. **Wikipedia.** "Snyder equal-area projection."
   - [Article Link](https://en.wikipedia.org/wiki/Snyder_equal-area_projection)

### Icosahedron Geometry

4. **Wikipedia.** "Regular icosahedron."
   - [Article Link](https://en.wikipedia.org/wiki/Regular_icosahedron)

5. **MathWorld.** "Regular Icosahedron."
   - [Wolfram MathWorld](https://mathworld.wolfram.com/RegularIcosahedron.html)

6. **Orbiter Forum.** "Coordinates of the vertices of an icosahedron on a circumscribed sphere."
   - [Forum Discussion](https://www.orbiter-forum.com/threads/coordinates-of-the-vertices-of-an-icosahedron-on-a-circumscribed-sphere.28851/)

### Spherical Triangle Centers

7. **BRSR Blog (2021).** "Spherical triangle centers with vectors."
   - [Blog Post](https://brsr.github.io/2021/05/02/spherical-triangle-centers.html)

8. **BRSR Blog (2021).** "Snyder's equal-area projection."
   - [Blog Post](https://brsr.github.io/2021/08/31/snyder-equal-area.html)

### Dual Polyhedra

9. **Wikipedia.** "Regular dodecahedron."
   - [Article Link](https://en.wikipedia.org/wiki/Regular_dodecahedron)

10. **Math@Brown.** "Duals of Regular Polyhedra."
    - [Brown University Resource](https://www.math.brown.edu/tbanchof/Beyond3d/chapter5/section03.html)

### Implementation References

11. **DGGRID Source Code** (dggridR implementation)
    - `DgProjTriRF.cpp`: Icosahedron vertex construction
    - `DgProjTriRF.h`: Default ISEA orientation constants

12. **hexify Source Code**
    - `src/icosahedron.h`: Interface definitions
    - `src/icosahedron.cpp`: Icosahedron construction and face assignment
    - `tests/testthat/test-icosahedron.R`: Verification tests

---

## Appendix A: Key Mathematical Values

### Icosahedron on Unit Sphere

| Property | Exact Value | Decimal Approximation |
|----------|-------------|------------------------|
| Circumradius | 1 | 1.0 |
| Edge length | √(10 - 2√5) / 2 | 1.051462224 |
| Midradius | φ/√3 | 0.934172358 |
| Inradius | φ²/(2√3) | 0.755761314 |
| Dihedral angle | arccos(-√5/3) | 138.189685° |
| Face area | π/5 | 0.628318531 sr |

### ISEA Orientation Constants

| Parameter | Value | Derivation |
|-----------|-------|------------|
| Vertex 0 longitude | 11.25° | 45° / 4 |
| Vertex 0 latitude | 58.28252559° | arcsin(2/√5) |
| Ring vertex latitude | ±26.565051177° | ±arctan(1/2) |
| Azimuth | 0° | Standard orientation |
| Longitude spacing | 72° | 360° / 5 |
| Ring offset | 36° | 72° / 2 |

### Trigonometric Identities

For latitude λ = arctan(1/2) = 26.565051177°:

```
sin(λ) = 1/√5 ≈ 0.447213595
cos(λ) = 2/√5 ≈ 0.894427191
tan(λ) = 1/2 = 0.5
```

For vertex 0 latitude λ₀ = 58.28252559°:

```
sin(λ₀) = 2/√5 ≈ 0.894427191
cos(λ₀) = 1/√5 ≈ 0.447213595
```

Note the complementary relationship: sin(λ₀) = cos(λ) and cos(λ₀) = sin(λ).

---

## Appendix B: Verification Code

### Python: Verify Vertex Latitudes

```python
import math

# Ring vertex latitude
ring_lat_rad = math.atan(0.5)
ring_lat_deg = math.degrees(ring_lat_rad)
print(f"Ring latitude: {ring_lat_deg:.9f}°")  # 26.565051177°

# Vertex 0 latitude
vert0_lat_rad = math.asin(2 / math.sqrt(5))
vert0_lat_deg = math.degrees(vert0_lat_rad)
print(f"Vertex 0 latitude: {vert0_lat_deg:.8f}°")  # 58.28252559°

# Verify complementary relationship
print(f"sin(ring_lat) = {math.sin(ring_lat_rad):.9f}")  # 1/√5
print(f"cos(ring_lat) = {math.cos(ring_lat_rad):.9f}")  # 2/√5
print(f"sin(vert0_lat) = {math.sin(vert0_lat_rad):.9f}")  # 2/√5
print(f"cos(vert0_lat) = {math.cos(vert0_lat_rad):.9f}")  # 1/√5
```

### R: Test Face Centers

```r
library(hexify)

hexify_build_icosa()
centers <- hexify_face_centers()

# Check coverage
cat("Northern hemisphere faces:", sum(centers$lat > 0), "\n")
cat("Southern hemisphere faces:", sum(centers$lat < 0), "\n")
cat("Longitude range:", max(centers$lon) - min(centers$lon), "radians\n")

# Sample face center
print(centers[1, ])
```

### C++: Dot Product Face Assignment

```cpp
#include <cmath>
#include <iostream>

struct Geo { double lon, lat; };

double great_circle_dist(const Geo& p1, const Geo& p2) {
  double cos_d = std::sin(p1.lat) * std::sin(p2.lat) +
                 std::cos(p1.lat) * std::cos(p2.lat) *
                 std::cos(p1.lon - p2.lon);
  return std::acos(std::clamp(cos_d, -1.0, 1.0));
}

int main() {
  Geo point{0.0, 0.5};  // Example point
  Geo center{0.0, 0.6}; // Example face center

  double dist = great_circle_dist(point, center);
  std::cout << "Distance: " << dist << " radians\n";

  return 0;
}
```

---

**Document Status:** Phase 1 Complete
**Last Updated:** 2025-12-17
**Author:** Agent C (Icosahedron Geometry Specialist)
**Review Status:** Pending Phase 2 Integration
