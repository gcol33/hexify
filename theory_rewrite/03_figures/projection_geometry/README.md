# Projection Geometry Figures

Mathematical illustrations for Lambert azimuthal equal-area projection and ISEA icosahedral face projection.

## Quick Start

Generate all figures:

```r
source("projection_diagrams.R")
generate_all_projection_diagrams(".")
```

## Files

- **projection_diagrams.R**: Complete R code for all diagrams
- **figure_spec.md**: Detailed specifications for each figure
- **fig_lambert_chord_distance.svg/.png**: Lambert chord distance geometry
- **fig_icosa_face_projection.svg/.png**: Icosahedron face projection concept

## Figures

### Figure 1: Lambert Chord Distance
Cross-section view showing the relationship between angular distance on a sphere and chord distance on a tangent plane. Formula: d = 2R sin(φ/2).

### Figure 2: Icosahedron Face Projection
Simplified view of sphere-to-plane projection for a single icosahedral face.

## Color Palette

- **Gray30** (`gray(0.3)`): Main structural elements
- **Gray50** (`gray(0.5)`): Secondary elements
- **Red** (`#E63946`): Points on sphere, distances
- **Blue** (`#457B9D`): Projected points, projection rays

## Requirements

Base R graphics only (no external packages required).
