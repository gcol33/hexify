## Lambert Azimuthal Equal-Area Projection

### Definition

The Lambert azimuthal equal-area projection maps a sphere of radius $R$ to a tangent plane while preserving area. For a projection centered at point $S$ with spherical coordinates $(\lambda_0, \phi_1)$, any point $P$ at coordinates $(\lambda, \phi)$ is mapped to planar coordinates $(x, y)$ such that infinitesimal areas on the sphere equal their corresponding areas in the plane (Snyder, 1987, p. 182).

The projection can be characterized geometrically through chord distance: a point $P$ at angular distance $c$ from the center projects to distance $\rho = 2R\sin(c/2)$ from the origin on the plane, measured along the azimuthal direction from $S$ to $P$ (Snyder, 1987, p. 182). This chord distance $\rho$ equals the straight-line distance through 3D space from $S$ to $P$ on the sphere's surface, derived from the law of cosines in the isosceles triangle formed by the sphere's center and the two surface points.

The projection is **not** perspective-based; it is a synthetic mathematical construction designed specifically to achieve equal-area and azimuthal properties simultaneously (Snyder, 1987, p. 182).

### Properties

**Equal-area preservation:** The projection maintains area exactly at all locations. This is the defining property distinguishing it from other azimuthal projections (Snyder, 1987, p. 188).

**Azimuthal from center:** All directions (azimuths) measured from the projection center are true. Straight lines radiating from the center represent great circle paths at their initial azimuth (Snyder, 1987, p. 182).

**Distortion characteristics:** The projection distorts shapes and angles increasingly with distance from center. At angular distance $c$ from center, the radial scale factor is $h' = \cos(c/2)$ (compression) and the tangential scale factor is $k' = \sec(c/2)$ (expansion), with $h' \cdot k' = 1$ enforcing area preservation (Snyder, 1987, eq. 24-22, 24-23, p. 188).

**Domain limitation:** The projection is bijective between the sphere minus the antipode point (angular distance $180°$ from center) and the open disk of radius $2R$ in the plane. The antipode maps to infinity and is typically excluded from the projection domain.

### Forward Formulas

For an oblique aspect projection centered at $(\lambda_0, \phi_1)$, the forward mapping $(\lambda, \phi) \to (x, y)$ is given by (Snyder, 1987, eq. 24-2 to 24-4, p. 185):

$$k' = \sqrt{\frac{2}{1 + \sin\phi_1\sin\phi + \cos\phi_1\cos\phi\cos(\lambda - \lambda_0)}}$$

$$x = R \cdot k' \cdot \cos\phi \cdot \sin(\lambda - \lambda_0)$$

$$y = R \cdot k' \cdot [\cos\phi_1\sin\phi - \sin\phi_1\cos\phi\cos(\lambda - \lambda_0)]$$

The factor $k'$ represents the local scale adjustment ensuring equal-area property while maintaining azimuthal directions.

For the **polar aspect** with center at the North Pole ($\phi_1 = 90°$), formulas simplify to (Snyder, 1987, eq. 24-5, 24-6, p. 185):

$$x = 2R\sin\left(\frac{\pi}{4} - \frac{\phi}{2}\right)\sin(\lambda - \lambda_0)$$

$$y = -2R\sin\left(\frac{\pi}{4} - \frac{\phi}{2}\right)\cos(\lambda - \lambda_0)$$

Here, $\rho = 2R\sin(\pi/4 - \phi/2) = 2R\sin((\pi/2 - \phi)/2)$ is the distance from origin, where $(\pi/2 - \phi)$ is the colatitude (angular distance from the North Pole).

### Inverse Formulas

The inverse mapping $(x, y) \to (\lambda, \phi)$ for the oblique aspect proceeds by first computing (Snyder, 1987, eq. 24-14 to 24-16, p. 187):

$$\rho = \sqrt{x^2 + y^2}$$

$$c = 2\arcsin\left(\frac{\rho}{2R}\right)$$

where $c$ is the angular distance from the projection center. Then:

$$\phi = \arcsin\left[\cos c \cdot \sin\phi_1 + \frac{y \cdot \sin c \cdot \cos\phi_1}{\rho}\right]$$

$$\lambda = \lambda_0 + \arctan\left[\frac{x \cdot \sin c}{\rho\cos\phi_1\cos c - y\sin\phi_1\sin c}\right]$$

If $\rho = 0$, the point is the projection center: $\phi = \phi_1, \lambda = \lambda_0$.

For the **polar aspect** (North Pole), the inverse simplifies to (Snyder, 1987, eq. 24-17 to 24-20, p. 187):

$$\phi = \frac{\pi}{2} - c = \frac{\pi}{2} - 2\arcsin\left(\frac{\rho}{2R}\right)$$

$$\lambda = \lambda_0 + \arctan\left(\frac{x}{-y}\right)$$

### Equal-Area Proof Sketch

**Scale factor approach:** For any map projection, let $h'$ be the scale factor in the direction of meridians (radial from center for azimuthal projections) and $k'$ the scale factor perpendicular to meridians. A projection is equal-area if and only if $h' \cdot k' = 1$ at every point (Snyder, 1987, p. 20-21).

For the Lambert azimuthal projection, at angular distance $c$ from center (Snyder, 1987, eq. 24-22, 24-23, p. 188):

$$h' = \cos\left(\frac{c}{2}\right), \quad k' = \sec\left(\frac{c}{2}\right) = \frac{1}{\cos(c/2)}$$

Therefore:
$$h' \cdot k' = \cos\left(\frac{c}{2}\right) \cdot \frac{1}{\cos(c/2)} = 1$$

This confirms the equal-area property at all finite points.

**Jacobian determinant verification (polar aspect):** For a coordinate transformation $(\lambda, \phi) \to (x, y)$, the area element transforms as:

$$dA_\text{plane} = \left|\det\left(\frac{\partial(x,y)}{\partial(\lambda,\phi)}\right)\right| d\lambda\, d\phi$$

On a sphere of radius $R$, the spherical area element is $dA_\text{sphere} = R^2\cos\phi\, d\lambda\, d\phi$. Equal-area projection requires $|\det(J)| = R^2\cos\phi$.

For the polar aspect, using forward formulas $x = 2R\sin(\pi/4 - \phi/2)\sin(\lambda - \lambda_0)$ and $y = -2R\sin(\pi/4 - \phi/2)\cos(\lambda - \lambda_0)$:

The Jacobian determinant evaluates to (calculation via chain rule and product rule):

$$\det(J) = \frac{\partial x}{\partial\lambda} \cdot \frac{\partial y}{\partial\phi} - \frac{\partial x}{\partial\phi} \cdot \frac{\partial y}{\partial\lambda}$$

Substituting partial derivatives and applying trigonometric identities yields:

$$\det(J) = 2R^2\sin\left(\frac{\pi}{4} - \frac{\phi}{2}\right)\cos\left(\frac{\pi}{4} - \frac{\phi}{2}\right)$$

Using $2\sin\alpha\cos\alpha = \sin(2\alpha)$:

$$\det(J) = R^2\sin\left(\frac{\pi}{2} - \phi\right) = R^2\cos\phi$$

This confirms the equal-area condition for the polar aspect. The oblique aspect follows similar principles with more complex algebra.

### Limitations

**Antipode singularity:** Points at angular distance $180°$ from the projection center (the antipodal point) cannot be projected; they map to infinity. Practical implementations restrict the projection domain to angular distances less than $\pi$.

**Increasing distortion:** Shape distortion increases with distance from center. The maximum angular distortion $\omega$ at angular distance $c$ is given by (Snyder, 1987, eq. 24-24, p. 188):

$$\sin\left(\frac{\omega}{2}\right) = \frac{k'^2 - 1}{k'^2 + 1}$$

At $c = 90°$ (quarter sphere), $k' = \sqrt{2}$ and $\omega \approx 70.5°$, indicating substantial shape distortion at hemisphere edges.

**Hemisphere constraint:** While the projection is mathematically defined for the entire sphere minus the antipode, practical applications typically restrict to a hemisphere or less to maintain acceptable shape fidelity. ISEA discrete global grids use this projection on icosahedral faces subtending approximately $72°$ from face center, where distortion remains moderate.

**No conformal property:** Unlike the stereographic projection, the Lambert azimuthal does not preserve angles. No projection can be both equal-area and conformal simultaneously—a fundamental constraint from differential geometry, as the sphere has positive Gaussian curvature while the plane has zero curvature (Snyder, 1987, p. 16-18).

### Conventions

**Coordinate ranges:**
- Latitude $\phi$: $-\pi/2 \le \phi \le \pi/2$ (radians), positive north
- Longitude $\lambda$: $-\pi \le \lambda \le \pi$ (radians), positive east

**Projection center:** $(\lambda_0, \phi_1)$ where $\lambda_0$ is the central meridian and $\phi_1$ is the latitude of projection center.

**Sphere radius:** Denoted $R$. For Earth applications, the authalic radius (radius of equal-area sphere matching Earth's surface area) is approximately $R_q = 6{,}371{,}007$ m for the WGS84 ellipsoid (Snyder, 1987, p. 187-188).

**Angular distance convention:** $c$ denotes angular distance along the great circle from projection center to point, with $0 \le c < \pi$.

### Sanity Check: Forward-Inverse Composition

The following R code verifies that composing forward and inverse projections returns the original coordinates (within numerical precision):

```r
# Forward projection (polar aspect, R=1, North Pole center)
forward_polar <- function(lon, lat) {
  rho <- 2 * sin(pi/4 - lat/2)
  x <- rho * sin(lon)
  y <- -rho * cos(lon)
  c(x = x, y = y)
}

# Inverse projection (polar aspect, R=1)
inverse_polar <- function(x, y) {
  rho <- sqrt(x^2 + y^2)
  c_dist <- 2 * asin(rho / 2)  # angular distance from pole
  lat <- pi/2 - c_dist
  lon <- atan2(x, -y)
  c(lon = lon, lat = lat)
}

# Test point: 45°E, 60°N
lon_in <- pi/4
lat_in <- pi/3
cat(sprintf("Input: lon = %.6f, lat = %.6f\n", lon_in, lat_in))

# Forward then inverse
xy <- forward_polar(lon_in, lat_in)
cat(sprintf("Projected: x = %.6f, y = %.6f\n", xy["x"], xy["y"]))

lonlat_out <- inverse_polar(xy["x"], xy["y"])
cat(sprintf("Recovered: lon = %.6f, lat = %.6f\n",
            lonlat_out["lon"], lonlat_out["lat"]))

# Verify round-trip error
error <- sqrt((lon_in - lonlat_out["lon"])^2 + (lat_in - lonlat_out["lat"])^2)
cat(sprintf("Round-trip error: %.2e radians\n", error))
stopifnot(error < 1e-12)
```

Expected output shows round-trip error below $10^{-12}$ radians, confirming numerical stability of the projection formulas.

---

**References:**

Snyder, J. P. (1987). *Map Projections: A Working Manual*. U.S. Geological Survey Professional Paper 1395. https://pubs.usgs.gov/pp/1395/report.pdf
