# Projection Geometry Figure Specifications

This document describes the mathematical figures for projection geometry in the hexify theory documentation.

## Color Palette

| Color | Hex Code | Usage |
|-------|----------|-------|
| Gray30 | `gray(0.3)` | Main structural elements (sphere, faces, axes) |
| Gray50 | `gray(0.5)` | Secondary elements (radii, grid lines, annotations) |
| Red | `#E63946` | Point P on sphere, chord distance |
| Blue | `#457B9D` | Projected point P', projection rays |

## Figure 1: Lambert Chord Distance Diagram

**Filename:** `fig_lambert_chord_distance.svg` / `.png`

**Purpose:** Illustrate the geometric relationship between angular distance on a sphere and chord distance on a tangent plane in Lambert azimuthal equal-area projection.

**Mathematical Content:**

- **Sphere:** Unit circle (R = 1) in cross-section, centered at origin O
- **Tangent plane:** Horizontal line at y = R, tangent to sphere at point S
- **Point S:** Tangent point at (0, R), top of sphere
- **Point P:** Point on sphere at angular distance φ from S
  - Coordinates: (R sin φ, R cos φ)
  - Default: φ = π/4 (45 degrees)
- **Point P':** Projection of P onto tangent plane
  - Lies on ray from O through P
  - Coordinates: (R sin φ / cos φ, R)
- **Chord distance d:** Distance from S to P' along tangent plane
  - Formula: d = 2R sin(φ/2)
  - Shown as thick red line

**Visual Elements:**

1. **Main structure** (Gray30):
   - Sphere outline (circle)
   - Tangent plane (horizontal line at y = R)
   - Label "Tangent plane" at left

2. **Construction lines** (Gray50, dashed):
   - Radius OS (vertical from center to tangent point)
   - Radius OP (from center to point P)

3. **Projection** (Blue):
   - Arrow from O through P to P'
   - Shows direction of projection

4. **Chord distance** (Red):
   - Thick line from S to P'
   - Label "d" at midpoint

5. **Angular distance** (Gray50):
   - Arc on sphere from S to P
   - Small arc near center showing angle φ
   - Label "φ" near angular arc

6. **Points:**
   - O (Gray30): Sphere center
   - S (Gray30): Tangent point
   - P (Red): Point on sphere
   - P' (Blue): Projected point on plane

7. **Formula annotation:**
   - d = 2R sin(φ/2) displayed below diagram

**Aspect ratio:** 1:1 (square)

**Dimensions:** 8 × 8 inches (SVG), 800 × 800 pixels (PNG)

---

## Figure 2: Icosahedron Face Projection

**Filename:** `fig_icosa_face_projection.svg` / `.png`

**Purpose:** Illustrate the concept of projecting points from the sphere surface onto an icosahedral face (tangent plane).

**Mathematical Content:**

- **Icosahedron face:** Equilateral triangle representing one of 20 icosahedral faces
  - Vertices V1, V2, V3
  - Centered at origin
  - Base horizontal
  - Light gray fill (90% opacity) to distinguish plane from background
- **Face center:** Centroid of triangle, marked with point
- **Sphere arc:** Curved arc above the triangle representing portion of sphere
  - Arc radius > triangle size to show curvature
  - Arc spans approximately the width of the triangle
- **Projection arrows:** Blue arrows from sphere to plane
  - Primary arrow (solid, thick): Main example showing point P → P'
  - Secondary arrow (dashed, thinner): Additional example to reinforce concept

**Visual Elements:**

1. **Tangent plane** (Gray30 outline, light gray fill):
   - Equilateral triangle with side length 3 units
   - Height = (√3/2) × 3 ≈ 2.6 units
   - Vertices labeled V1, V2, V3 (Gray50)
   - Annotation: "Tangent plane (icosa face)"

2. **Face center** (Gray30):
   - Point at triangle centroid
   - Label "Face center"

3. **Sphere surface** (Gray30):
   - Arc curving above triangle
   - Thick line (3 pt)
   - Annotation: "Sphere surface"

4. **Projection 1** (primary example):
   - Point P (Red): Point on sphere arc, left of center
   - Point P' (Blue): Projected point on face
   - Arrow (Blue, solid, 2.5 pt): From P to P'

5. **Projection 2** (secondary example):
   - Unnamed points on sphere and face
   - Arrow (Blue, dashed, 1.5 pt): From sphere to face
   - Shows projection is not just one point

6. **Scale indicator** (Gray50):
   - Vertical line at left showing distance R
   - Marked with 0 at bottom, R at top
   - Label "Scale" rotated 90°

**Aspect ratio:** 1:1 (square)

**Dimensions:** 8 × 8 inches (SVG), 800 × 800 pixels (PNG)

**Plot window:** xlim = (-2, 2), ylim = (-1.5, 2.5)

---

## Usage Instructions

### Generate all figures

```r
source("projection_diagrams.R")
generate_all_projection_diagrams("output_directory")
```

### Generate individual figures

```r
source("projection_diagrams.R")

# Lambert chord distance (default φ = π/4)
create_lambert_chord_diagram("path/to/fig_lambert_chord_distance")

# Lambert chord distance with custom angle
create_lambert_chord_diagram("path/to/fig_lambert_chord_distance", phi = pi/3)

# Icosahedron face projection
create_icosa_face_projection_diagram("path/to/fig_icosa_face_projection")
```

### Preview without saving

```r
source("projection_diagrams.R")

# Preview in R graphics device
create_lambert_chord_diagram()  # output_path = NULL
create_icosa_face_projection_diagram()  # output_path = NULL
```

---

## Design Rationale

### Visual Hierarchy

1. **Primary focus** (Red): Points and distances being measured/projected
2. **Projection mechanism** (Blue): Arrows and projected results
3. **Structure** (Gray30): Main geometric objects (sphere, faces)
4. **Construction** (Gray50): Supporting elements (radii, arcs, labels)

### Simplifications

- Cross-section view for Lambert diagram reduces 3D → 2D while preserving key geometry
- Single icosahedral face shown (not full icosahedron) to avoid visual clutter
- Scale indicator instead of axis labels for cleaner presentation

### Mathematical Accuracy

- Lambert chord distance formula d = 2R sin(φ/2) is exact
- Sphere arc and tangent plane relationship is geometrically correct
- Projection rays follow true gnomonic projection (through sphere center)
- All angles and distances are computed, not hand-drawn

---

## Output Files

When `generate_all_projection_diagrams()` is run, the following files are created:

```
output_directory/
├── fig_lambert_chord_distance.svg
├── fig_lambert_chord_distance.png
├── fig_icosa_face_projection.svg
└── fig_icosa_face_projection.png
```

**SVG files:** Vector graphics, suitable for LaTeX/PDF documents, scalable without quality loss

**PNG files:** Raster graphics at 100 DPI (800×800 px), suitable for web/presentations

---

## Quality Assurance Checklist

- [x] Geometrically accurate representations (formulas verified)
- [x] Clear visual hierarchy (color-coded by function)
- [x] Minimal color palette (4 colors: 2 grays + red + blue)
- [x] Consistent styling across both figures
- [x] Labels positioned to avoid overlaps
- [x] High-resolution output (800×800 px for PNG)
- [x] Vector output available (SVG)
- [x] Base R graphics (no external dependencies)
- [x] Reproducible code with clear helper functions
