## Icosahedron Geometry

The Icosahedral Snyder Equal Area (ISEA) projection begins by mapping the sphere to a regular icosahedron, which is then unfolded into planar triangular faces. This section establishes the geometric foundations: what an icosahedron is, how it is oriented, how face centers are computed, and how points are assigned to faces.

### Definition and Properties

A regular icosahedron is a convex Platonic solid with 12 vertices, 30 edges of equal length, and 20 equilateral triangular faces. When inscribed in a unit sphere (radius $R = 1$), each spherical triangle face covers an area of $\pi/5 \approx 0.628$ steradians, ensuring the 20 faces partition the sphere uniformly.

The icosahedron exhibits the highest symmetry of any polyhedron compatible with spherical partitioning: icosahedral symmetry (group $I_h$). Vertex coordinates can be expressed using the golden ratio $\phi = (1 + \sqrt{5})/2 \approx 1.618$. The 12 vertices consist of all cyclic permutations of $(0, \pm 1, \pm \phi)$, normalized to lie on the unit sphere.

| Property | Value | Derivation |
|----------|-------|------------|
| Vertices | 12 | Fundamental to icosahedron |
| Edges | 30 | Euler characteristic: $V - E + F = 2$ |
| Faces | 20 | Equal-area spherical triangles |
| Face area | $\pi/5$ | $4\pi / 20$ steradians |
| Inradius | $\phi^2/(2\sqrt{3}) \approx 0.756$ | Distance to face center |
| Edge length | $\sqrt{10 - 2\sqrt{5}}/2 \approx 1.051$ | On unit sphere |

### Standard ISEA Orientation

The standard ISEA orientation, established by Snyder (1992) and adopted by PROJ, DGGRID, dggridR, and hexify, places the icosahedron with:

**Vertex 0:**
$$\text{lon} = 11.25^\circ \text{ E}, \quad \text{lat} = 58.28252559^\circ \text{ N}$$

**Azimuth:** $0^\circ$ (no rotation around the polar axis)

The remaining 11 vertices form two pentagonal rings plus a south pole vertex. The longitude $11.25^\circ = 45^\circ/4$ was chosen for computational convenience and to minimize distortion over major landmasses.

#### Vertex 0 Latitude Derivation

The latitude $58.28252559^\circ$ arises from the requirement that the icosahedron be properly oriented with respect to the sphere. Specifically:

$$\text{lat}_0 = \arcsin\left(\frac{2}{\sqrt{5}}\right) \approx 58.28252559^\circ$$

This can also be expressed as $90^\circ - \arccos(1/\sqrt{5})$. The exact trigonometric values are:

$$\sin(\text{lat}_0) = \frac{2}{\sqrt{5}} \approx 0.8944, \quad \cos(\text{lat}_0) = \frac{1}{\sqrt{5}} \approx 0.4472$$

These values ensure the icosahedron's vertices align correctly with the sphere's curvature, maintaining equal edge lengths and proper face symmetry.

### Pentagon Locations and Ring Vertices

When one vertex is placed at the north pole, the remaining 10 middle vertices arrange into two pentagonal rings. These rings occur at latitudes:

$$\text{lat}_{\text{ring}} = \pm \arctan\left(\frac{1}{2}\right) = \pm 26.565051177^\circ$$

#### Arctan(1/2) Derivation

This value emerges from the icosahedron's intrinsic geometry. When a vertex occupies the north pole, the ratio of vertical to horizontal distance for the adjacent ring vertices is exactly 1:2. This gives:

$$\tan(\text{colatitude}) = 2 \quad \Rightarrow \quad \text{colatitude} = \arctan(2) \approx 63.435^\circ$$

Therefore:
$$\text{latitude} = 90^\circ - \arctan(2) = \arctan\left(\frac{1}{2}\right) \approx 26.565^\circ$$

Equivalently:
$$\sin(\text{lat}_{\text{ring}}) = \frac{1}{\sqrt{5}} \approx 0.4472, \quad \cos(\text{lat}_{\text{ring}}) = \frac{2}{\sqrt{5}} \approx 0.8944$$

Note the complementary relationship with Vertex 0: $\sin(\text{lat}_0) = \cos(\text{lat}_{\text{ring}})$ and vice versa.

**Ring structure:**
- Upper ring (5 vertices): $+26.565^\circ$ latitude, longitudes spaced at $72^\circ$ intervals
- Lower ring (5 vertices): $-26.565^\circ$ latitude, offset by $36^\circ$ from upper ring
- South pole vertex: $-90^\circ$ latitude

The $72^\circ$ spacing derives from $360^\circ/5$, and the $36^\circ$ offset ($72^\circ/2$) creates the pentagonal gyroelongated bipyramid structure characteristic of the icosahedron.

| Parameter | Value | Derivation |
|-----------|-------|------------|
| Vertex 0 longitude | $11.25^\circ$ | $45^\circ / 4$ |
| Vertex 0 latitude | $58.28252559^\circ$ | $\arcsin(2/\sqrt{5})$ |
| Ring latitudes | $\pm 26.565051177^\circ$ | $\pm \arctan(1/2)$ |
| Pentagon spacing | $72^\circ$ | $360^\circ / 5$ |
| Ring offset | $36^\circ$ | $72^\circ / 2$ |

### Face Centers

Each icosahedral face is a spherical triangle. The face center is defined as the centroid of this spherical triangle, computed using the Euclidean centroid method:

1. Convert the three vertices $v_1, v_2, v_3$ to Cartesian coordinates on the unit sphere
2. Sum the vertex vectors: $\mathbf{v}_{\text{sum}} = v_1 + v_2 + v_3$
3. Normalize to unit length: $\mathbf{v}_{\text{center}} = \mathbf{v}_{\text{sum}} / |\mathbf{v}_{\text{sum}}|$
4. Convert back to geographic coordinates: $\text{lon} = \arctan_2(y, x)$, $\text{lat} = \arcsin(z)$

This centroid represents the point of maximum symmetry within each face. For equilateral spherical triangles, this coincides with other classical triangle centers (circumcenter, incenter). Unlike planar geometry, the spherical centroid does not generally divide the triangle into three equal-area sub-triangles.

The face centers are precomputed once during initialization. For the standard ISEA orientation, the 20 face centers span both hemispheres and distribute evenly across longitude, ensuring global coverage.

### Face Assignment Algorithm

**Problem:** Given a point $(\text{lon}, \text{lat})$ on the sphere, determine which of the 20 icosahedral faces contains it.

**Solution:** Use the great-circle distance from the point to each face center. The face with the minimum distance (equivalently, maximum dot product) contains the point.

#### Mathematical Justification

The great-circle distance $d$ between two points on a unit sphere is given by the spherical law of cosines:

$$\cos(d) = \sin(\text{lat}_1) \cdot \sin(\text{lat}_2) + \cos(\text{lat}_1) \cdot \cos(\text{lat}_2) \cdot \cos(\text{lon}_1 - \text{lon}_2)$$

Converting both points to Cartesian coordinates $\mathbf{p} = (x_1, y_1, z_1)$ and $\mathbf{c} = (x_2, y_2, z_2)$, the dot product is:

$$\mathbf{p} \cdot \mathbf{c} = x_1 x_2 + y_1 y_2 + z_1 z_2 = \cos(d)$$

Since cosine is monotonically decreasing on $[0, \pi]$, maximizing $\mathbf{p} \cdot \mathbf{c}$ minimizes $d$. Each face occupies a region bounded by three great circle arcs, and the face center is the most representative point of that region. For points strictly inside a face (not on boundaries), the nearest face center is always the center of the containing face.

#### Implementation

The algorithm iterates through all 20 face centers, computing the great-circle distance to each:

```cpp
int best = 0;
double best_dist = great_circle_distance(point, center[0]);

for (int i = 1; i < 20; ++i) {
  double dist = great_circle_distance(point, center[i]);
  if (dist < best_dist) {
    best = i;
    best_dist = dist;
  }
}
return best;
```

An optimization is possible: since cosine is monotonic, comparing $\cos(d)$ directly (larger is better) avoids the `acos()` call, improving performance.

### Edge Cases

**Points on face boundaries:** When a point lies exactly on an edge shared by two faces, both face centers are equidistant. The algorithm deterministically returns the first face encountered in iteration order (typically the lowest face number). Since edges have measure zero (probability 0 of being hit by a random point), this ambiguity has negligible practical impact.

**Points on vertices:** Three faces share each vertex. The algorithm assigns the point to one of them (the first in iteration order). Again, this is deterministic and consistent within each implementation.

**Numerical precision:** Floating-point arithmetic introduces small errors. Values that are mathematically equal may differ slightly in computation. The implementation clamps intermediate values to $[-1, 1]$ to ensure trigonometric functions remain well-defined.

**Uniqueness guarantee:** For points strictly inside a face, the assignment is unique and mathematically well-defined. On boundaries, the assignment is ambiguous geometrically but deterministic algorithmically. Different implementations may assign boundary points differently, but the same implementation is always consistent. Since the icosahedron serves as an intermediate structure (points are further subdivided into hexagonal cells), boundary ambiguity does not affect final cell assignments in practice.

### References

Snyder, J.P. (1992). An Equal-Area Map Projection For Polyhedral Globes. *Cartographica*, 29(1):10-21.

PROJ Contributors. Icosahedral Snyder Equal Area. PROJ documentation. https://proj.org/en/stable/operations/projections/isea.html

Coxeter, H.S.M. (1973). *Regular Polytopes*. Dover Publications. Chapter 10: The Regular Polyhedra.
