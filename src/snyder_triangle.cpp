#include "snyder_triangle.h"
#include <cmath>
#include <algorithm>

namespace hexify {

namespace {

inline double dot(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void cross(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

inline double norm(const double a[3]) { return std::sqrt(dot(a, a)); }

inline double det3(const double a[3], const double b[3], const double c[3]) {
  double bc[3];
  cross(b, c, bc);
  return dot(a, bc);
}

// The arc between two unit vectors, accurate at every angle.
inline double arc(const double a[3], const double b[3]) {
  double ab[3];
  cross(a, b, ab);
  return std::atan2(norm(ab), dot(a, b));
}

// Signed area of the spherical triangle (a, b, c), positive counterclockwise
// seen from outside: tan(A / 2) = det / (1 + a.b + b.c + c.a).
inline double signed_area(const double a[3], const double b[3], const double c[3]) {
  return 2.0 * std::atan2(det3(a, b, c), 1.0 + dot(a, b) + dot(b, c) + dot(c, a));
}

} // anon

SnyderTriangle make_snyder_triangle(const double v0[3], const double v1[3],
                                    const double v2[3]) {
  SnyderTriangle t;
  for (int k = 0; k < 3; ++k) {
    t.v0[k] = v0[k];
    t.v1[k] = v1[k];
    t.v2[k] = v2[k];
  }
  t.area = signed_area(v0, v1, v2);
  t.det = det3(v0, v1, v2);
  t.c01 = dot(v0, v1);
  t.c12 = dot(v1, v2);
  t.c20 = dot(v2, v0);
  double x12[3];
  cross(v1, v2, x12);
  t.s12 = norm(x12);
  t.theta12 = std::atan2(t.s12, t.c12);
  return t;
}

void snyder_triangle_forward(const SnyderTriangle& t, const double v[3],
                             double& b1, double& b2) {
  // p: where the great circle from v0 through v meets the edge v1 v2.
  double n0[3], n12[3], p[3];
  cross(t.v0, v, n0);
  cross(t.v1, t.v2, n12);
  cross(n0, n12, p);
  const double pn = norm(p);
  if (pn == 0.0) {  // v = v0
    b1 = b2 = 0.0;
    return;
  }
  const double sgn = (dot(p, t.v1) + dot(p, t.v2) < 0.0) ? -1.0 : 1.0;
  for (int k = 0; k < 3; ++k) p[k] *= sgn / pn;

  // h^2 = (1 - v0.v) / (1 - v0.p), as a ratio of half-arc sines
  const double h = std::sin(0.5 * arc(t.v0, v)) / std::sin(0.5 * arc(t.v0, p));
  b2 = h * signed_area(t.v0, t.v1, p) / t.area;
  b1 = h - b2;
}

double snyder_triangle_inverse_ray(const SnyderTriangle& t, double b1, double b2,
                                   double p[3]) {
  const double h = b1 + b2;

  // a = A(v0, v1, p): the area fraction b2 / h of the triangle. The arc from
  // v1 to p follows from tan of half of it = g / f.
  const double a = (h > 0.0) ? (b2 / h) * t.area : 0.0;
  const double S = std::sin(a);
  const double sh = std::sin(0.5 * a);
  const double C = 2.0 * sh * sh;  // 1 - cos(a)
  const double f = S * t.det + C * (t.c01 * t.c12 - t.c20);
  const double g = C * t.s12 * (1.0 + t.c01);
  const double s1p = 2.0 * std::atan2(g, f);

  // p = Slerp(v1, v2), at arc s1p from v1
  const double w1 = std::sin(t.theta12 - s1p) / t.s12;
  const double w2 = std::sin(s1p) / t.s12;
  for (int k = 0; k < 3; ++k) p[k] = w1 * t.v1[k] + w2 * t.v2[k];

  // sin(z / 2) = h sin(z_p / 2), z_p the arc v0 p
  return 2.0 * std::asin(std::min(1.0, h * std::sin(0.5 * arc(t.v0, p))));
}

void snyder_triangle_inverse(const SnyderTriangle& t, double b1, double b2,
                             double out[3]) {
  double p[3];
  const double z = snyder_triangle_inverse_ray(t, b1, b2, p);
  const double cp = dot(t.v0, p);
  double u[3];
  for (int k = 0; k < 3; ++k) u[k] = p[k] - cp * t.v0[k];
  const double un = norm(u);
  const double cz = std::cos(z), sz = std::sin(z);
  for (int k = 0; k < 3; ++k) out[k] = cz * t.v0[k] + sz * u[k] / un;
}

} // namespace hexify
