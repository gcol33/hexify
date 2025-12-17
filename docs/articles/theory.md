# Mathematical Foundations

## Overview

hexify implements the **ISEA (Icosahedral Snyder Equal Area)** discrete
global grid system. This vignette explains the mathematical foundations:
the projection geometry, coordinate systems, aperture subdivision, and
cell indexing.

## The Problem: Equal-Area Grids on a Sphere

Any projection from a sphere to a plane must distort *something*. For
spatial statistics, we need cells of equal area regardless of location.
Standard latitude-longitude grids fail badly: a 1° cell at the equator
covers ~12,300 km², while the same cell near the poles covers a tiny
fraction of that.

The ISEA projection solves this by: 1. Inscribing a regular icosahedron
in the sphere 2. Projecting each spherical “cap” onto its corresponding
flat triangular face using a modified Lambert equal-area projection 3.
Overlaying a hexagonal grid on the resulting planar triangles

## The Lambert Azimuthal Equal-Area Projection

The foundation of Snyder’s projection is Lambert’s azimuthal equal-area
projection, developed by Johann Heinrich Lambert in 1772.

### Geometric Definition

The key insight is elegantly simple. Consider a sphere with a plane
tangent at point S:

![](theory_files/figure-html/lambert-geometry-1.svg)

**The Lambert projection rule:** For any point P on the sphere, measure
the straight-line distance d through 3D space from the tangent point S
to P (the chord distance). Place the projected point P’ on the plane at
distance d from S, in the same azimuthal direction.

### Why Chord Distance Preserves Area

For a point at angular distance $`\phi`$ from S (measured along the
sphere surface), the chord distance is:

``` math
d = 2R \sin\left(\frac{\phi}{2}\right)
```

This specific relationship ensures that a small area element $`dA`$ on
the sphere maps to an equal area $`dA`$ on the plane. The
area-preserving property arises because the Jacobian determinant of this
transformation equals 1.

![](theory_files/figure-html/lambert-area-preservation-1.svg)

Each colored band has equal area on the sphere. After Lambert
projection, the shapes change (outer bands stretch radially, compress
tangentially) but the areas remain equal. This is the defining property
of an equal-area projection.

**Reference:** [Wolfram MathWorld: Lambert Azimuthal Equal-Area
Projection](https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html)

## From Lambert to Snyder: The Icosahedron

Lambert’s projection works for a single tangent plane covering at most a
hemisphere. To cover the entire globe with minimal distortion, Snyder’s
key insight was to use **20 tangent planes**—one for each face of a
regular icosahedron.

### The Icosahedron

A regular icosahedron has: - **20 triangular faces** (each covering
~1/20 of Earth’s surface ≈ 25.5 million km²) - **12 vertices** (where 5
faces meet—these become pentagon cells) - **30 edges**

![](theory_files/figure-html/icosahedron-concept-1.svg)

Because each face subtends a relatively small solid angle (~1/20 of the
sphere), projection distortion within each face is limited. The maximum
angular distortion is bounded at approximately 17.3°.

### Standard ISEA Orientation

The icosahedron orientation is chosen to place vertices (pentagon cells)
in oceanic or low-population areas:

    #> Linking to GEOS 3.13.1, GDAL 3.11.0, PROJ 9.6.0; sf_use_s2() is TRUE
    #> Warning in st_point_on_surface.sfc(sf::st_zm(x)): st_point_on_surface may not
    #> give correct results for longitude/latitude data

![](theory_files/figure-html/face-centers-1.svg)

**Reference:** [Sahr, White & Kimerling (2003). Geodesic discrete global
grid
systems.](https://www.researchgate.net/publication/246557072_ISEA_discrete_global_grids)

## Snyder’s Projection Formulas

Snyder adapted Lambert’s approach for triangular faces. The projection
involves:

1.  Finding which face contains the point
2.  Computing spherical coordinates relative to the face center
3.  Applying transformation formulas that preserve area

### The Forward Projection

For a point P at geographic coordinates (lon, lat), relative to a face
with center at (lon₀, lat₀):

**Step 1:** Compute angular distance z and azimuth Az from face center:
``` math
z = \arccos(\sin\phi_0 \sin\phi + \cos\phi_0 \cos\phi \cos(\lambda - \lambda_0))
```
``` math
Az = \arctan\left(\frac{\cos\phi \sin(\lambda - \lambda_0)}{\cos\phi_0 \sin\phi - \sin\phi_0 \cos\phi \cos(\lambda - \lambda_0)}\right)
```

**Step 2:** Reduce azimuth to \[0°, 120°) sector and apply Snyder’s
transformation:
``` math
\delta_z = \arctan\left(\frac{\tan E_l}{\cos Az + \cot 30° \sin Az}\right)
```
``` math
h = \arccos(\sin Az \sin G \cos E_l - \cos Az \cos G)
```
``` math
A_G = Az + G + h - \pi
```
``` math
Az' = \arctan\left(\frac{2 A_G}{R_1^2 \tan^2 E_l - 2 A_G \cot 30°}\right)
```

**Step 3:** Compute radial distance:
``` math
f = \frac{\tan E_l}{2(\cos Az' + \cot 30° \sin Az') \sin(\delta_z/2)}
```
``` math
\rho = 2 R_1 f \sin(z/2)
```

**Step 4:** Convert to Cartesian:
``` math
x = \rho \sin Az' + x_0, \quad y = \rho \cos Az' + y_0
```

### Key Constants

| Constant | Value        | Description                          |
|----------|--------------|--------------------------------------|
| $`E_l`$  | 37.37736814° | Face center to edge angular distance |
| G        | 36°          | Geometric angle (360°/10)            |
| $`R_1`$  | 0.9103832815 | Scaling factor for equal-area        |

**Reference:** [Snyder (1992). An equal-area map projection for
polyhedral
globes.](https://www.researchgate.net/publication/296900899_Snyder_equal-area_map_projection_for_polyhedral_globes)

### The Inverse Projection (Newton-Raphson Iteration)

The inverse—converting planar coordinates back to geographic—cannot be
solved analytically because the forward transformation involves
transcendental functions. Instead, a Newton-Raphson iterative method is
used.

![](theory_files/figure-html/newton-raphson-1.svg)

The iteration solves for the azimuth Az given the planar coordinates.
Starting from an initial estimate, each iteration refines:

``` math
Az_{n+1} = Az_n - \frac{f(Az_n)}{f'(Az_n)}
```

Convergence typically occurs within 3-5 iterations to machine precision.

**Reference:** [PROJ Documentation: Icosahedral Snyder Equal
Area](https://proj.org/en/stable/operations/projections/isea.html)

## Aperture and Resolution

**Aperture** determines how cells subdivide at each resolution
level—it’s the ratio of parent cell area to child cell area.

![](theory_files/figure-html/aperture-diagrams-1.svg)

| Aperture | Children/Parent | Area Ratio | Best For |
|----|----|----|----|
| 3 | 3 | 3:1 | Fine resolution control, dggridR compatibility |
| 4 | 4 | 4:1 | Power-of-2 scaling, GIS workflows |
| 7 | 7 | 7:1 | Rapid cell growth, coarse analysis |

### Cell Count Growth

``` r

cat("Resolution  Aperture 3    Aperture 4    Aperture 7\n")
#> Resolution  Aperture 3    Aperture 4    Aperture 7
cat("---------  ----------    ----------    ----------\n")
#> ---------  ----------    ----------    ----------
for (res in 0:8) {
  cells_ap3 <- 10 * 3^res + 2
  cells_ap4 <- 10 * 4^res + 2
  cells_ap7 <- 10 * 7^res + 2
  cat(sprintf("    %d      %10s    %10s    %10s\n",
              res,
              format(cells_ap3, big.mark = ","),
              format(cells_ap4, big.mark = ","),
              format(cells_ap7, big.mark = ",")))
}
#>     0              12            12            12
#>     1              32            42            72
#>     2              92           162           492
#>     3             272           642         3,432
#>     4             812         2,562        24,012
#>     5           2,432        10,242       168,072
#>     6           7,292        40,962     1,176,492
#>     7          21,872       163,842     8,235,432
#>     8          65,612       655,362    57,648,012
```

**Reference:** [Location coding on icosahedral aperture 3 hexagon
discrete global
grids](https://www.researchgate.net/publication/223763428_Location_coding_on_icosahedral_aperture_3_hexagon_discrete_global_grids)

## Orientation Classes

The aperture determines how child cells align with parent cells:

![](theory_files/figure-html/orientation-classes-1.svg)

For **aperture 3**, orientation alternates between Class I and II at
each resolution level. For **aperture 7**, each level adds a rotation of
arctan(√3/5) ≈ 19.1°.

## Pentagon Cells

At every resolution, exactly **12 cells are pentagons** (not hexagons).
These occur at the icosahedron vertices and have area exactly 5/6 that
of hexagonal cells.

![](theory_files/figure-html/pentagon-locations-1.svg)

| Location   | Latitude | Longitudes                  |
|------------|----------|-----------------------------|
| Poles      | ±90°     | 0°                          |
| Upper ring | 26.57°   | 0°, 72°, 144°, 216°, 288°   |
| Lower ring | -26.57°  | 36°, 108°, 180°, 252°, 324° |

The latitude 26.57° = arctan(1/2) arises from the geometry of a regular
icosahedron.

## Coordinate Systems

hexify uses several coordinate systems internally:

| System | Components | Description |
|----|----|----|
| **GEO** | lon, lat | WGS84 degrees |
| **Icosa Triangle** | face (0-19), tx, ty | Triangle face + normalized coords \[0,1\] |
| **Quad IJ** | quad (0-11), i, j | Quad + integer grid indices |
| **SEQNUM** | integer | Global cell ID (1-based, dggridR compatible) |

### Triangle to Quad

The 20 triangular faces are paired into 12 “quads” (rhombi). Each quad
combines two adjacent triangles sharing an edge, simplifying the grid
indexing.

## Round-Trip Accuracy

``` r

# Round-trip accuracy test
original_lon <- 16.37
original_lat <- 48.21

cat(sprintf("Original: (%.4f, %.4f)\n\n", original_lon, original_lat))
#> Original: (16.3700, 48.2100)

for (ap in c(3, 4, 7)) {
  res <- if (ap == 7) 6 else 10
  grid <- hex_grid(resolution = res, aperture = ap)
  cell_id <- lonlat_to_cell(original_lon, original_lat, grid)
  recovered <- cell_to_lonlat(cell_id, grid)
  error_km <- sqrt((recovered$lon - original_lon)^2 +
                   (recovered$lat - original_lat)^2) * 111
  cat(sprintf("Aperture %d (res %2d): cell %d -> (%.4f, %.4f), ~%.1f km from center\n",
              ap, res, cell_id,
              recovered$lon, recovered$lat, error_km))
}
#> Aperture 3 (res 10): cell 126594 -> (16.4204, 48.2877), ~10.3 km from center
#> Aperture 4 (res 10): cell 2245587 -> (16.4177, 48.2119), ~5.3 km from center
#> Aperture 7 (res  6): cell 237704 -> (14.7883, 46.5600), ~253.7 km from center
```

## Summary of ISEA Properties

| Property | Description |
|----|----|
| **Equal area** | All hexagonal cells identical; pentagons = 5/6 hex area |
| **Bounded distortion** | Max angular distortion ≈ 17.3° at face edges |
| **Uniform topology** | Hexagons have 6 neighbors; pentagons have 5 |
| **No perfect nesting** | Cells at different resolutions partially overlap |

## References

- Snyder, J. P. (1992). [An equal-area map projection for polyhedral
  globes](https://www.researchgate.net/publication/296900899_Snyder_equal-area_map_projection_for_polyhedral_globes).
  *Cartographica*, 29(1), 10-21.

- Sahr, K., White, D., & Kimerling, A. J. (2003). [Geodesic discrete
  global grid
  systems](https://www.researchgate.net/publication/246557072_ISEA_discrete_global_grids).
  *Cartography and Geographic Information Science*, 30(2), 121-134.

- [Wolfram MathWorld: Lambert Azimuthal Equal-Area
  Projection](https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html)

- [PROJ Documentation: Icosahedral Snyder Equal
  Area](https://proj.org/en/stable/operations/projections/isea.html)

- [DGGRID software](https://github.com/sahrk/DGGRID) by Kevin Sahr
