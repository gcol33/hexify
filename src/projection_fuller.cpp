#include "projection_fuller.h"
#include "constants.h"
#include <cmath>

namespace hexify {

namespace {

const double SQRT3 = std::sqrt(3.0);

// The plane triangle of Gray's equations (6)-(8) is centred on its centroid
// with edge ARC; the face plane has unit edge and its origin at the lower-left
// vertex.
const double ORIGIN_X = 0.5;
const double ORIGIN_Y = 1.0 / (2.0 * SQRT3);

} // anon

std::pair<double,double> fuller_face_xy(double z, double az) {
  // The point on the secondary plane triangle (Gray eqs. 14-16), with the
  // azimuth measured from the vertex SV1 on the +y axis.
  const double r = kFullerZ0 * std::tan(z);
  const double xs = r * std::sin(az);
  const double ys = r * std::cos(az);

  // Lengths along the secondary triangle's edges (Gray eqs. 20, 22, 24),
  // each less EL/2.
  const double a1s = 2.0 * ys / SQRT3 + kFullerEL / 3.0 - 0.5 * kFullerEL;
  const double a2s = xs - ys / SQRT3 + kFullerEL / 3.0 - 0.5 * kFullerEL;
  const double a3s = kFullerEL / 3.0 - xs - ys / SQRT3 - 0.5 * kFullerEL;

  // Arc lengths a_i - alpha along the spherical triangle's edges
  // (Gray eqs. 29-31).
  const double b1 = std::atan(a1s / kFullerDVE);
  const double b2 = std::atan(a2s / kFullerDVE);
  const double b3 = std::atan(a3s / kFullerDVE);

  // The plane point (Gray eqs. 36, 38).
  const double xpp = 0.5 * (b2 - b3);
  const double ypp = (2.0 * b1 - b2 - b3) / (2.0 * SQRT3);

  return {xpp / kFullerArc + ORIGIN_X, ypp / kFullerArc + ORIGIN_Y};
}

std::pair<double,double> fuller_face_polar(double x, double y, double tol,
                                           int max_iters, int* iters) {
  const double xpp = (x - ORIGIN_X) * kFullerArc;
  const double ypp = (y - ORIGIN_Y) * kFullerArc;

  // Gray's equation (39) in Crider's form (27): with u = a2 - alpha,
  //   tan(u - A) + tan(u) + tan(u - B) + tan(alpha) = 0,
  // A = x'' - sqrt(3) y'', B = 2 x''. Every term increases with u, so the
  // root is unique and Newton's method converges from the root of the
  // linearised equation.
  const double A = xpp - SQRT3 * ypp;
  const double B = 2.0 * xpp;
  double u = (A + B - kFullerTanAlpha) / 3.0;
  int k = 0;
  for (; k < max_iters; ++k) {
    const double t1 = std::tan(u - A);
    const double t2 = std::tan(u);
    const double t3 = std::tan(u - B);
    const double g = t1 + t2 + t3 + kFullerTanAlpha;
    const double dg = 3.0 + t1 * t1 + t2 * t2 + t3 * t3;
    const double du = g / dg;
    u -= du;
    if (std::fabs(du) <= tol) { ++k; break; }
  }
  if (iters) *iters = k;

  // a1 - alpha = u - A (Crider eq. 35); back to the secondary plane triangle
  // (Crider eqs. 36-38).
  const double a1s = kFullerDVE * std::tan(u - A);
  const double a2s = kFullerDVE * std::tan(u);
  const double ys = 0.5 * SQRT3 * (a1s + kFullerEL / 2.0 - kFullerEL / 3.0);
  const double xs = a2s + kFullerEL / 2.0 + ys / SQRT3 - kFullerEL / 3.0;

  // Onto the sphere along the ray through the plane point (Crider eq. 39).
  const double z = std::atan2(std::hypot(xs, ys), kFullerZ0);
  double az = std::atan2(xs, ys);
  if (az < 0.0) az += kTwoPi;
  return {z, az};
}

} // namespace hexify
