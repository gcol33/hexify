#include "projection_ivea.h"
#include "constants.h"
#include "dual.h"
#include <algorithm>
#include <cmath>

namespace hexify {

namespace {

// The unit-edge plane triangle about its centre C': the vertex B' of the
// canonical right triangle lies at the circumradius along +y, and the edge
// midpoint A' at the inradius along the unit vector (sqrt(3)/2, 1/2), azimuth
// 60 degrees. Azimuths are read atan2(x, y), as Snyder's.
constexpr double kRv = 0.57735026918962576451;    // 1 / sqrt(3)
constexpr double kRin = 0.28867513459481288225;   // 1 / (2 sqrt(3))
constexpr double kUx = 0.86602540378443864676;    // sqrt(3) / 2
constexpr double kUy = 0.5;

// cos and sin of k * 120 degrees, k = 0, 1, 2
constexpr double kSectorCos[3] = {1.0, -0.5, -0.5};
constexpr double kSectorSin[3] = {0.0, 0.86602540378443864676, -0.86602540378443864676};

} // anon

template <class T>
std::pair<T,T> ivea_face_xy(const VertexGcParams& p, T z, T az) {
  using std::acos;
  using std::atan;
  using std::atan2;
  using std::cos;
  using std::sin;
  using std::sqrt;

  // The right triangle holding the point: sector k runs from vertex k
  // (azimuth 0) to vertex k + 1 (azimuth 120 degrees) with its edge midpoint
  // at 60; past the midpoint the triangle is the mirror image about it.
  const int k = azimuth_sector(value_of(az));
  T a = az - k * k2PiOver3;
  const bool mirror = value_of(a) > kPiOver3;
  if (mirror) a = k2PiOver3 - a;

  // In the triangle C B P: CB = BC, CP = z, angle a at C. sin(x / 2), x = BP,
  // by the haversine formula, and the angle theta at B from BC to BP.
  const T sz = sin(z), cz = cos(z);
  const T hz = sin((z - p.bc) * 0.5);
  const T ha = sin(a * 0.5);
  const T sin_half_x = sqrt(hz * hz + p.sin_bc * sz * ha * ha);
  const T theta = atan2(sin(a) * sz, p.sin_bc * cz - p.cos_bc * sz * cos(a));

  // The cut through B and P at angle rho from BA meets CA at D, at angle
  // delta: cos(delta) = sin(rho) cos(AB) (eq. 22). The triangle B C D holds
  // the share (beta + gamma - rho - delta) / excess of the area (eq. 20),
  // which is C'D' / C'A' on the plane (eqs. 21, 23).
  const T rho = p.beta - theta;
  const T delta = acos(sin(rho) * p.cos_ab);
  const T frac = (p.beta + p.gamma - rho - delta) / p.excess;
  // Right angle at A: tan(BD) = tan(AB) / cos(rho)
  const T bd = atan(p.tan_ab / cos(rho));
  // B'P' / B'D' = sin(x / 2) / sin(BD / 2) (eq. 28)
  const T s = sin_half_x / sin(bd * 0.5);

  const T dx = frac * (kRin * kUx);
  const T dy = frac * (kRin * kUy);
  T qx = s * dx;
  T qy = kRv + s * (dy - kRv);
  if (mirror) {
    const T w = 2.0 * (qx * kUx + qy * kUy);
    qx = w * kUx - qx;
    qy = w * kUy - qy;
  }
  const double c = kSectorCos[k], sn = kSectorSin[k];
  return {qx * c + qy * sn + 0.5, qy * c - qx * sn + kRin};
}

template std::pair<double,double> ivea_face_xy<double>(const VertexGcParams&, double, double);
template std::pair<Dual2,Dual2> ivea_face_xy<Dual2>(const VertexGcParams&, Dual2, Dual2);

std::pair<double,double> ivea_face_polar(const VertexGcParams& p, double x, double y) {
  const double X = x - 0.5;
  const double Y = y - kRin;
  double phi = std::atan2(X, Y);
  if (phi < 0.0) phi += kTwoPi;
  const int k = azimuth_sector(phi);

  // Turned back into sector 0, and mirrored into the canonical triangle
  const double c = kSectorCos[k], sn = kSectorSin[k];
  double qx = X * c - Y * sn;
  double qy = Y * c + X * sn;
  const bool mirror = phi - k * k2PiOver3 > kPiOver3;
  if (mirror) {
    const double w = 2.0 * (qx * kUx + qy * kUy);
    qx = w * kUx - qx;
    qy = w * kUy - qy;
  }

  // The line from B' = (0, kRv) through the point meets C'A' at D' = t u.
  const double wx = qx, wy = qy - kRv;
  const double w_len = std::hypot(wx, wy);
  double z, a;
  if (w_len < 1e-15) {
    z = p.bc;
    a = 0.0;
  } else {
    // B' + s w = t u, crossed with w: t = (B' x w) / (u x w)
    const double t = (-kRv * wx) / (kUx * wy - kUy * wx);
    const double frac = t / kRin;
    const double s = w_len / std::hypot(t * kUx, t * kUy - kRv);

    // rho + delta = S, with cos(delta) = sin(rho) cos(AB):
    //   cos(S - rho) = sin(rho) cos(AB)  =>  tan(rho) = cos(S) / (cos(AB) - sin(S)),
    // and rho in [0, beta] fixes the branch.
    const double S = p.beta + p.gamma - frac * p.excess;
    const double rho = std::atan2(-std::cos(S), std::sin(S) - p.cos_ab);
    const double bd = std::atan(p.tan_ab / std::cos(rho));
    const double sin_half_x = s * std::sin(0.5 * bd);
    const double xb = 2.0 * std::asin(std::min(1.0, sin_half_x));
    const double theta = p.beta - rho;

    // In the triangle B C P: BC, BP = xb, angle theta at B
    const double sx = std::sin(xb), cx = std::cos(xb);
    const double hb = std::sin(0.5 * (p.bc - xb));
    const double ht = std::sin(0.5 * theta);
    z = 2.0 * std::asin(std::min(1.0, std::sqrt(hb * hb + p.sin_bc * sx * ht * ht)));
    a = std::atan2(std::sin(theta) * sx, p.sin_bc * cx - p.cos_bc * sx * std::cos(theta));
  }
  if (mirror) a = k2PiOver3 - a;
  return {z, a + k * k2PiOver3};
}

std::pair<double,double> ivea_face_polar_newton(const VertexGcParams& p, double x, double y) {
  // The step after one of size d is of size ~d^2, so a step below 1e-14
  // leaves the point at full double precision.
  constexpr double kTol = 1e-14;
  constexpr int kMaxIters = 100;

  // Start from the plane point's own polar form, with the circumradius
  // stretched to the arc BC.
  const double X = x - 0.5, Y = y - kRin;
  double z = std::hypot(X, Y) / kRv * p.bc;
  double az = std::atan2(X, Y);
  if (az < 0.0) az += kTwoPi;
  for (int i = 0; i < kMaxIters; ++i) {
    const auto f = ivea_face_xy(p, Dual2(z, 1.0, 0.0), Dual2(az, 0.0, 1.0));
    const double rx = f.first.v - x, ry = f.second.v - y;
    // At a vertex, where the start can land, the derivative is singular.
    if (std::hypot(rx, ry) < 1e-15) break;
    const double a11 = f.first.d[0], a12 = f.first.d[1];
    const double a21 = f.second.d[0], a22 = f.second.d[1];
    const double det = a11 * a22 - a12 * a21;
    const double dz = (a22 * rx - a12 * ry) / det;
    const double da = (a11 * ry - a21 * rx) / det;
    z -= dz;
    az -= da;
    if (az < 0.0) az += kTwoPi;
    if (az >= kTwoPi) az -= kTwoPi;
    if (std::fabs(dz) + std::fabs(da) <= kTol) break;
  }
  return {z, az};
}

} // namespace hexify
