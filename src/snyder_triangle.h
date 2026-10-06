#pragma once
// Snyder's (1992) equal-area map between one spherical triangle and one plane
// triangle, in the unit-vector form of Recht (2021). The spherical triangle
// (v0, v1, v2) maps to the plane triangle with barycentric coordinates
// (b0, b1, b2); lines from v0 map to lines from the plane image of v0, and the
// area on one side of such a line keeps its fraction of the triangle's area.
// Both directions are closed-form. The triangle may be irregular and either
// orientation.

namespace hexify {

struct SnyderTriangle {
  double v0[3], v1[3], v2[3];
  double area;     // signed spherical area A(v0, v1, v2)
  double det;      // det[v0, v1, v2]
  double c01, c12, c20;
  double s12;      // |v1 x v2|
  double theta12;  // arc v1 v2
};

SnyderTriangle make_snyder_triangle(const double v0[3], const double v1[3],
                                    const double v2[3]);

// Unit vector v -> barycentric (b1, b2); b0 = 1 - b1 - b2.
void snyder_triangle_forward(const SnyderTriangle& t, const double v[3],
                             double& b1, double& b2);

// Barycentric (b1, b2) -> unit vector.
void snyder_triangle_inverse(const SnyderTriangle& t, double b1, double b2,
                             double out[3]);

// The inverse as a ray from v0: returns the arc z from v0 to the point and
// sets p to where its great circle from v0 meets the edge v1 v2. The point
// lies at arc z from v0 towards p.
double snyder_triangle_inverse_ray(const SnyderTriangle& t, double b1, double b2,
                                   double p[3]);

} // namespace hexify
