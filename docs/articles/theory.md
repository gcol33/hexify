# Mathematical Foundations

## Overview

hexify implements the **ISEA (Icosahedral Snyder Equal Area)** discrete
global grid system. This vignette explains the mathematical foundations
with intuitive geometric illustrations.

## The Problem: Mapping a Sphere to a Plane

The fundamental challenge in cartography is mapping the curved surface
of Earth onto a flat plane. Any such mapping must distort
*something*—you cannot perfectly preserve both area and shape
simultaneously. For ecological and statistical applications,
**equal-area** projections are essential: we need grid cells of equal
size regardless of location.

### Why Equal-Area Matters

Consider counting species in grid cells:

- If cells vary in area, counts aren’t comparable
- A cell twice as large will, on average, contain twice as many
  individuals
- Standard lat/lon grids fail badly—cells near poles are tiny compared
  to equatorial cells

The ISEA projection solves this by preserving area while minimizing
shape distortion.

## The Lambert Azimuthal Equal-Area Projection

Before understanding Snyder’s projection, we need to understand its
foundation: the **Lambert azimuthal equal-area projection**, developed
by Johann Heinrich Lambert in 1772.

### The Core Idea: Chord Distance

The key insight is beautifully simple. Imagine a sphere with a flat
plane tangent to it at point S (the “south pole” in polar aspect):

![Lambert projection: Points on the sphere are mapped to the plane using
chord distance.](theory_files/figure-html/lambert-concept-1.svg)

Lambert projection: Points on the sphere are mapped to the plane using
chord distance.

**How it works:**

1.  Pick any point **P** on the sphere
2.  Measure the straight-line (chord) distance **d** from **S** to **P**
3.  Place the projected point **P’** on the plane at distance **d** from
    **S**, in the same direction

This simple rule—*use chord distance*—automatically preserves area! The
mathematical magic is that the area element transforms correctly.

### Why Does This Preserve Area?

The chord distance formula for a point at angular distance φ from S is:

``` math
d = 2R \sin\left(\frac{\phi}{2}\right)
```

where R is the sphere’s radius. This specific relationship ensures that
a small patch of area dA on the sphere maps to an equal area dA on the
plane.

![Equal-area property: Concentric circles on the sphere map to circles
on the plane with areas
preserved.](theory_files/figure-html/lambert-area-1.svg)

Equal-area property: Concentric circles on the sphere map to circles on
the plane with areas preserved.

The colored bands have equal area on the sphere. After Lambert
projection, the shapes change (bands near the edge are stretched
radially, compressed tangentially), but the *areas* remain equal.

## From Lambert to Snyder: Projecting onto an Icosahedron

Lambert’s projection works for a single tangent plane. But what if we
want to cover the entire globe with minimal distortion? Snyder’s
insight: use **20 tangent planes**, one for each face of an icosahedron.

### The Icosahedron Advantage

![An icosahedron inscribed in a sphere provides 20 nearly-flat
triangular faces.](theory_files/figure-html/icosahedron-concept-1.svg)

An icosahedron inscribed in a sphere provides 20 nearly-flat triangular
faces.

Why an icosahedron?

1.  **20 identical faces**: Each triangular face covers ~1/20 of Earth
    (about 25.5 million km²)
2.  **Low distortion**: Because each face is small and nearly flat,
    projection distortion stays below 17.3°
3.  **Triangles subdivide naturally**: Each triangle can be recursively
    divided into smaller triangles

### Snyder’s Projection: Modified Lambert for Triangles

Snyder adapted Lambert’s method for projecting spherical triangles (the
curved patches on the globe) to planar triangles (the flat icosahedron
faces):

![Snyder's projection maps curved spherical triangles to flat planar
triangles while preserving
area.](theory_files/figure-html/snyder-process-1.svg)

Snyder’s projection maps curved spherical triangles to flat planar
triangles while preserving area.

**The Snyder projection process:**

1.  **Subdivide each spherical face** into 6 smaller triangles (from
    center to vertices)
2.  **Apply modified Lambert** to each sub-triangle, preserving area
3.  **Reassemble** the planar triangles to form the flat icosahedron
    face

The mathematical details involve computing intersection points of great
circles and carefully chosen scaling factors, but the principle is the
same as Lambert: preserve area by using the right distance relationship.

## The Complete Projection Pipeline

Now let’s see how hexify implements this:

### Step 1: Determine the Face

For any point on Earth, first find which of the 20 icosahedral faces
contains it. The algorithm computes the angular distance from the point
to each face center and selects the closest face.

### Step 2: Project to Face Coordinates

Apply Snyder’s equal-area formulas to compute (x, y) coordinates on the
planar triangle. The formulas involve:

- Computing the azimuth angle from face center to point
- Computing the angular distance from face center
- Applying the area-preserving transformation

### Step 3: Normalize and Quantize

Convert the (x, y) coordinates to normalized (tx, ty) in range \[0, 1\],
then quantize to grid cell indices based on the resolution.

``` r

library(hexify)

# Project a point through the pipeline
result <- hexify_forward(lon = 16.37, lat = 48.21)
cat(sprintf("Face: %d, tx: %.4f, ty: %.4f\n",
            result["face"], result["tx"], result["ty"]))
#> Face: 2, tx: NA, ty: NA
```

### Icosahedron Orientation

The icosahedron can be oriented arbitrarily on the globe. The default
ISEA orientation places vertex 0 at:

- Longitude: 45°/4
- Latitude: arctan(√5/2)

This orientation was chosen by Sahr et al. (2003) to minimize distortion
over land masses while placing the 12 pentagonal cells (which occur at
icosahedron vertices) in oceanic or low-population areas. The specific
vertex latitude arctan(√5/2) arises from the geometry of a regular
icosahedron inscribed in a sphere.

``` r

# Get all face centers
centers <- hexify_face_centers()
head(centers, 5)
#>          lon       lat
#> 1 -1.3744468 1.2059325
#> 2 -0.5890486 0.6154797
#> 3  0.1963495 0.3648638
#> 4  0.9817477 0.6154797
#> 5  1.7671459 1.2059325
```

## Hexagonal Grid Generation

### From Triangles to Hexagons

Each triangular face is subdivided into hexagonal cells:

![Hexagonal grid overlaid on a triangular face. Hexagons tile the plane
efficiently with equal-area
cells.](theory_files/figure-html/hexagon-subdivision-1.svg)

Hexagonal grid overlaid on a triangular face. Hexagons tile the plane
efficiently with equal-area cells.

The hexagonal grid provides several advantages:

1.  **Uniform neighbors**: Each hexagon has exactly 6 neighbors (unlike
    squares with 8)
2.  **Isotropic**: Distance to all neighbors is equal
3.  **Efficient packing**: Hexagons cover the plane with minimal
    perimeter per area

### Aperture and Resolution

**Aperture** determines how many child cells fit in one parent cell.
This is the key parameter controlling grid resolution:

![Aperture controls how parent cells subdivide into children. Aperture 3
divides into 3 children, aperture 4 into 4, and aperture 7 into
7.](theory_files/figure-html/aperture-visual-1.svg)

Aperture controls how parent cells subdivide into children. Aperture 3
divides into 3 children, aperture 4 into 4, and aperture 7 into 7.

| Aperture | Children per Parent | Class Pattern | Description |
|----|----|----|----|
| 3 | 3 | I/II alternating | Fine resolution control |
| 4 | 4 | I only | Power-of-2 scaling |
| 7 | 7 | III only | Fastest cell growth |
| 4/3 | 4 then 3 | Mixed | Aperture 4 for first m levels, then 3 |

**Resolution** is the number of subdivision levels. Cell count formulas:

| Aperture | Cell Count Formula       |
|----------|--------------------------|
| 3        | 10 × 3^res + 2           |
| 4        | 10 × 4^res + 2           |
| 7        | 10 × 7^res + 2           |
| 4/3      | 10 × 4^m × 3^(res-m) + 2 |

``` r

# Cell count growth for all apertures
cat("Resolution  Aperture 3    Aperture 4    Aperture 7\n")
#> Resolution  Aperture 3    Aperture 4    Aperture 7
cat("---------  ----------    ----------    ----------\n")
#> ---------  ----------    ----------    ----------
for (res in 0:5) {
  cells_ap3 <- 10 * 3^res + 2
  cells_ap4 <- 10 * 4^res + 2
  cells_ap7 <- 10 * 7^res + 2
  cat(sprintf("    %d      %9d     %9d     %9d\n",
              res, cells_ap3, cells_ap4, cells_ap7))
}
#>     0             12            12            12
#>     1             32            42            72
#>     2             92           162           492
#>     3            272           642          3432
#>     4            812          2562         24012
#>     5           2432         10242        168072
```

For **mixed aperture (4/3)**, the first `m` levels use aperture 4, then
remaining levels use aperture 3:

``` r

# Mixed aperture example: m=2 means first 2 levels use aperture 4
m <- 2
cat("Mixed 4/3 (m=2): first 2 levels aperture 4, rest aperture 3\n")
#> Mixed 4/3 (m=2): first 2 levels aperture 4, rest aperture 3
for (res in 0:5) {
  if (res <= m) {
    cells <- 10 * 4^res + 2
  } else {
    cells <- 10 * 4^m * 3^(res - m) + 2
  }
  cat(sprintf("  Res %d: %d cells\n", res, cells))
}
#>   Res 0: 12 cells
#>   Res 1: 42 cells
#>   Res 2: 162 cells
#>   Res 3: 482 cells
#>   Res 4: 1442 cells
#>   Res 5: 4322 cells
```

### Orientation Classes

DGGRID defines three orientation classes based on how child cells align
with parent cells during subdivision:

- **Class I**: Child cell vertices align with parent cell vertices.
  Hexagons have a “flat top” orientation.
- **Class II**: Child cell centers align with parent cell vertices.
  Hexagons are rotated 30° (“pointy top”).
- **Class III**: Neither alignment—child cells are rotated relative to
  parent cells (skewed grids).

Each aperture produces a characteristic class pattern:

| Aperture | Class Pattern | Description |
|----|----|----|
| 3 | I → II → I → II… | Alternates between Class I and II at each resolution |
| 4 | I → I → I… | Always Class I (aligned subdivision) |
| 7 | III → III → III… | Always Class III (arctan(√3/5) rotation per level) |
| 4/3 | I → I → … → I/II | Class I for aperture-4 levels, then alternates |

For **aperture 3**, odd resolutions (1, 3, 5, …) produce Class I grids
with flat-top hexagons, while even resolutions (2, 4, 6, …) produce
Class II grids with pointy-top hexagons.

For **aperture 4**, the 2×2 subdivision preserves alignment, so all
resolutions are Class I.

For **aperture 7**, the subdivision involves a √7 scaling and
arctan(√3/5) rotation, creating the characteristic Class III skewed
pattern. This rotation accumulates: resolution 2 is rotated
2×arctan(√3/5) from resolution 0, and so on.

For **mixed aperture (4/3)**, the class pattern follows aperture 4
(Class I) for the first `m` levels, then switches to aperture 3’s
alternating I/II pattern for subsequent levels.

The class affects hexagon orientation and the coordinate system
internals, but not the fundamental cell assignment or equal-area
property.

## Space-Filling Curves

### Why Space-Filling Curves?

Space-filling curves provide:

1.  **Compact representation**: Cell identity as a single integer/string
2.  **Hierarchical encoding**: Parent-child relationships from index
    structure
3.  **Spatial locality**: Nearby cells often have similar indices

Each aperture requires a different curve that matches its subdivision
geometry.

### Z3 Encoding (Aperture 3)

For aperture 3, the Z3 curve uses base-3 (ternary) encoding. Each
resolution level adds one ternary digit representing which of the 3
child cells contains the point.

    Index = FF D₁D₂...Dₙ
    where:
      FF = face number (00-19)
      Dᵢ = digit 0, 1, or 2 (base 3)

The three child positions form a triangular arrangement within each
parent cell. The digit values correspond to:

- **0**: Center-bottom child
- **1**: Upper-left child
- **2**: Upper-right child

This encoding preserves hierarchy: truncating the last digit gives the
parent cell index.

``` r

# Z3 index structure
idx <- hexify_lonlat_to_index(16.37, 48.21, resolution = 5, aperture = 3)
cat("Index:", idx, "\n")
#> Index: 0211101
cat("Face:", substr(idx, 1, 2), "\n")
#> Face: 02
cat("Path:", substr(idx, 3, nchar(idx)), "\n")
#> Path: 11101

# Parent is obtained by dropping the last digit
parent_idx <- substr(idx, 1, nchar(idx) - 1)
cat("Parent:", parent_idx, "\n")
#> Parent: 021110
```

### Z-Order (Morton) Encoding (Aperture 4)

For aperture 4, the Morton curve (Z-order curve) interleaves the bits of
i and j grid coordinates. This classic space-filling curve maps a 2D
grid to a 1D index while preserving locality.

The 2×2 subdivision at each level creates 4 child cells arranged in a
square:

      +---+---+
      | 2 | 3 |
      +---+---+
      | 0 | 1 |
      +---+---+

Each child position is encoded as 2 bits: the low bit is from i, the
high bit is from j. For resolution r, the index uses 2r bits total.

    For cell at grid position (i, j):
      bit 0 = i₀, bit 1 = j₀   (resolution 1)
      bit 2 = i₁, bit 3 = j₁   (resolution 2)
      ...

This interleaving creates the characteristic “Z” pattern that gives the
curve its name.

``` r

# Z-order index
idx4 <- hexify_lonlat_to_index(16.37, 48.21, resolution = 5, aperture = 4)
cat("Z-order index:", idx4, "\n")
#> Z-order index: 0233232

# The index encodes (quad, i, j) implicitly
coords <- hexify_index_to_cell(idx4, aperture = 4)
cat("Decoded - Quad:", coords$quad, "i:", coords$i, "j:", coords$j, "\n")
#> Decoded - Quad: i: 31 j: 26
```

### Z7 Encoding (Aperture 7)

Aperture 7 uses base-7 encoding with a twist: each subdivision involves
both scaling by √7 and rotation by arctan(√3/5). The 7 child cells form
a hexagonal rosette pattern:

          1
        6   2
          0
        5   3
          4

The center cell (0) is surrounded by 6 peripheral cells (1-6) arranged
in a hexagonal ring. Unlike apertures 3 and 4, the child cells are
rotated relative to the parent, which is why aperture 7 produces Class
III grids.

The encoding uses Eisenstein integers (complex numbers of the form a +
bω where ω = e^(2πi/3)) to handle the hexagonal geometry. Each digit 0-6
represents a position in the rosette.

``` r

# Z7 index
idx7 <- hexify_lonlat_to_index(16.37, 48.21, resolution = 3, aperture = 7)
cat("Z7 index:", idx7, "\n")
#> Z7 index: 02153

# Hierarchy still works - parent is obtained by dropping last digit
parent7 <- substr(idx7, 1, nchar(idx7) - 1)
cat("Parent:", parent7, "\n")
#> Parent: 0215
```

### Mixed Aperture (4/3) Encoding

Mixed aperture grids combine the aperture 4 and aperture 3 encoding
schemes. The first `m` resolution levels use Morton (Z-order) encoding,
then subsequent levels use Z3 encoding.

This hybrid approach provides:

- Faster initial cell growth (aperture 4 at coarse scales)
- Finer resolution control (aperture 3 at fine scales)
- dggridR compatibility with ISEA43H grids

``` r

# Mixed aperture index (via hexify() function)
result <- hexify(data.frame(lon = 16.37, lat = 48.21),
                 lon = "lon", lat = "lat", area = 1000, aperture = "4/3")
cat("Mixed 4/3 cell_id:", result$cell_id, "\n")
```

### Comparison of Encodings

| Property | Z3 (Aperture 3) | Morton (Aperture 4) | Z7 (Aperture 7) | Mixed (4/3) |
|----|----|----|----|----|
| Base | 3 (ternary) | 4 (quaternary) | 7 | 4 then 3 |
| Digits per level | 1 | 2 bits | 1 | Mixed |
| Child arrangement | Triangle | Square | Hex rosette | Square → Triangle |
| Rotation per level | 0° or 30° | 0° | arctan(√3/5) | Mixed |
| Locality | Good | Excellent | Good | Good |

## Inverse Projection

Converting cell coordinates back to lon/lat requires:

1.  **Decode index**: Extract face and (i, j) coordinates
2.  **Compute cell center**: Find center point in face coordinates
3.  **Inverse Snyder**: Convert face (x, y) to geographic (lon, lat)

``` r

# Round-trip test for all apertures
original_lon <- 16.37
original_lat <- 48.21

cat(sprintf("Original coordinates: (%.4f, %.4f)\n\n", original_lon, original_lat))
#> Original coordinates: (16.3700, 48.2100)

for (ap in c(3, 4, 7)) {
  # Use appropriate resolution for each aperture
  res <- if (ap == 7) 6 else 10
  idx <- hexify_lonlat_to_index(original_lon, original_lat,
                                 resolution = res, aperture = ap)
  recovered <- hexify_index_to_lonlat(idx, aperture = ap)
  error <- sqrt((recovered["lon"] - original_lon)^2 +
                (recovered["lat"] - original_lat)^2)
  cat(sprintf("Aperture %d (res %2d): recovered (%.4f, %.4f), error: %.6f°\n",
              ap, res, recovered["lon"], recovered["lat"], error))
}
#> Aperture 3 (res 10): recovered (16.4204, 48.2877), error: 0.092576°
#> Aperture 4 (res 10): recovered (16.4177, 48.2119), error: 0.047712°
#> Aperture 7 (res  6): recovered (-10.0330, 3.8584), error: 51.615693°
```

## Cell Properties

### Equal Area

A key property of ISEA grids: all cells have approximately equal area.

``` r

# Compare cell areas across apertures at similar resolutions
earth_area_km2 <- 510072000

cat("Aperture comparison at ~1000 km² target:\n\n")
#> Aperture comparison at ~1000 km² target:
for (ap in c(3, 4, 7)) {
  grid <- hexify_grid(area = 1000, aperture = ap)
  cat(sprintf("Aperture %d: resolution %d → %.1f km² (%.0f cells)\n",
              ap, grid$resolution, grid$area, grid$total_cells))
}

# Mixed aperture available via hexify() function, not hexify_grid()
cat("\n(Mixed aperture 4/3 is available via hexify() with aperture = '4/3')\n")
#> 
#> (Mixed aperture 4/3 is available via hexify() with aperture = '4/3')
```

### Cell Shape

Cells are approximately hexagonal, but:

- Cells at face edges may be irregular
- Cell shape varies slightly with latitude due to projection
- Pentagon cells exist at icosahedron vertices (12 total)

## References

- Sahr, K., White, D., & Kimerling, A. J. (2003). Geodesic discrete
  global grid systems.
- Snyder, J. P. (1992). An equal-area map projection for polyhedral
  globes.
- DGGRID documentation: <https://github.com/sahrk/DGGRID>
