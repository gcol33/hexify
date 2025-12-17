# Newton-Raphson Iteration Figure Specification

## Purpose

Visualizes the geometric interpretation of Newton-Raphson iteration for solving the Snyder azimuth residual equation in the inverse projection.

## Mathematical Context

The inverse projection requires solving the transcendental equation:

$$f(\text{Az}) = \text{agh} - \text{Az} - G + (\pi - h) = 0$$

where Az is the spherical azimuth. Newton-Raphson iteration starts with an initial guess (the planar azimuth) and converges to the root via:

$$\text{Az}_{n+1} = \text{Az}_n - \frac{f(\text{Az}_n)}{f'(\text{Az}_n)}$$

## Visual Elements

### Function Representation
- **Curve**: Representative function $f(x) = x^3 - 2x - 2$
  - Chosen to mimic Snyder residual properties: smooth, monotonic, single root
  - Color: Gray30 (#4D4D4D)

### Iteration Sequence
- **Initial guess** $x_0 = 2.5$: Starting point offset from root for visibility
- **Iterations** $x_1, x_2, x_3, x_4$: Successive approximations
- **Tangent lines**: Red dashed lines showing Newton step geometry
  - From $(x_n, f(x_n))$ to $(x_{n+1}, 0)$
  - Illustrates how derivative determines next guess

### Convergence Markers
- **Iteration points**: Red filled circles on curve
- **Converged point** $x_4$: Red cross (×) on x-axis
- **True root**: Black cross (×) on x-axis (reference)

### Axes and Grid
- **Horizontal line** at $y = 0$: Root-finding target
- **Vertical grid**: Light gray dotted lines for readability
- **Axis labels**: "Az (degrees)" (x), "f(Az)" (y)

## Color Palette

| Element | Color | Hex Code | Rationale |
|---------|-------|----------|-----------|
| Function curve | Gray30 | #4D4D4D | High contrast without distraction |
| Axes/grid | Gray50 | #808080 | Subtle reference structure |
| Iteration elements | Red | #E63946 | Accent for active process |
| True root | Black | #000000 | Ground truth reference |

## Convergence Demonstration

The figure illustrates **quadratic convergence**: each iteration approximately doubles the number of correct digits. Expected progression:

| Iteration | Error (approx) | Correct Digits |
|-----------|----------------|----------------|
| 0 | 0.7 | 0.1 |
| 1 | 0.05 | 1.3 |
| 2 | 0.002 | 2.7 |
| 3 | 3e-6 | 5.5 |
| 4 | 9e-12 | 11.0 |

This mirrors the empirical performance described in section text: 3-5 iterations to reach $10^{-12}$ tolerance.

## Usage Context

This figure supports the "Newton-Raphson Formulation" and "Convergence Properties" subsections of the Inverse Projection section. It provides visual intuition for:

1. Why the method works (tangent line approximation)
2. How quickly it converges (rapid reduction in error)
3. The geometric meaning of the iteration formula

## Output Files

- `newton_raphson.svg`: Vector format for publications
- `newton_raphson.png`: Raster format (300 DPI) for web/presentations

## Execution

Run from package root:
```r
source("theory_rewrite/03_figures/newton/newton_raphson.R")
```

Outputs generated in `theory_rewrite/03_figures/newton/`.

## Design Rationale

**Base R graphics**: Chosen for precise control over mathematical plotting elements and reproducibility.

**No rainbow colors**: Grayscale palette with single accent color maintains professional appearance and ensures accessibility (colorblind-safe).

**Iteration count**: Four iterations shown to demonstrate convergence without visual clutter. Matches typical performance (median: 4 iterations).

**Representative function**: $f(x) = x^3 - 2x - 2$ avoids exposing implementation details of the Snyder residual while preserving key mathematical properties (smooth, monotonic, transcendental).
