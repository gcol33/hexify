// projection_ak.cpp - Kaseorg's octahedral projection (see projection_ak.h)
//
// Written from Griffin (2026), Sec. 11a, and the forward formula of libhex9's
// core/h9_math.h (h9_ak_core), used as a reference only.

#include "projection_ak.h"
#include "constants.h"
#include "dual.h"
#include "hex9_warp.h"
#include <cmath>
#include <type_traits>

namespace hexify {

namespace {

// The face's vertices in its own frame: the centre on +z, the first vertex
// at azimuth 0 (+y), the second at 240 degrees and the third at 120, each at
// the arc g from the centre, cos g = 1/sqrt(3). On the octahedron they are
// orthonormal.
const double kInvSqrt2 = 0.70710678118654752440;
const double kInvSqrt3 = 0.57735026918962576451;
const double kInvSqrt6 = 0.40824829046386301637;
const double kVertex[3][3] = {
  {0.0, 2.0 * kInvSqrt6, kInvSqrt3},
  {-kInvSqrt2, -kInvSqrt6, kInvSqrt3},
  {kInvSqrt2, -kInvSqrt6, kInvSqrt3}
};

// Face-plane (x, y) as weights on the face's vertices: top, lower left,
// lower right
template <class T>
inline void plane_weights(T x, T y, T b[3]) {
  b[0] = y / kSin60;
  b[2] = x - 0.5 * b[0];
  b[1] = 1.0 - b[0] - b[2];
}

// The unnormalised direction of the point with weights b, along the vertices
template <class T>
inline void kaseorg(const T b[3], T X[3]) {
  using std::tan;
  using std::sqrt;
  T t[3], t2[3];
  for (int k = 0; k < 3; k++) {
    t[k] = tan((0.5 * kPi) * b[k]);
    t2[k] = t[k] * t[k];
  }
  for (int k = 0; k < 3; k++) {
    const T& a = t2[(k + 1) % 3];
    const T& c = t2[(k + 2) % 3];
    X[k] = t[k] * sqrt(sqrt(a + c + kKaseorgAlpha * a * c));
  }
}

// Face-plane (x, y) -> (z, az), az in (-pi, pi]
template <class T>
std::pair<T,T> polar_of_plane(T x, T y) {
  using std::atan2;
  using std::sqrt;
  T b[3], X[3];
  plane_weights(x, y, b);
  kaseorg(b, X);
  T p[3];
  for (int r = 0; r < 3; r++) {
    p[r] = X[0] * kVertex[0][r] + X[1] * kVertex[1][r] + X[2] * kVertex[2][r];
  }
  return {atan2(sqrt(p[0] * p[0] + p[1] * p[1]), p[2]), atan2(p[0], p[1])};
}

// (z, az) -> face-plane (x, y) by Newton's method on two weights. The point's
// components c_k along the vertices fix its direction; with m the largest,
// the direction of X(b) is c's where X_i c_m = X_m c_i for the other two.
std::pair<double,double> plane_of_polar(double z, double az) {
  const double p[3] = {std::sin(z) * std::sin(az), std::sin(z) * std::cos(az), std::cos(z)};
  double c[3];
  int m = 0;
  for (int k = 0; k < 3; k++) {
    c[k] = p[0] * kVertex[k][0] + p[1] * kVertex[k][1] + p[2] * kVertex[k][2];
    if (c[k] > c[m]) m = k;
  }
  // A vertex is its own corner.
  if (c[m] > 1.0 - 1e-15) {
    const double corner[3][2] = {{0.5, kSin60}, {0.0, 0.0}, {1.0, 0.0}};
    return {corner[m][0], corner[m][1]};
  }
  const double sum = c[0] + c[1] + c[2];
  double b0 = c[0] / sum, b2 = c[2] / sum;
  const int i1 = (m + 1) % 3, i2 = (m + 2) % 3;
  for (int iter = 0; iter < 50; iter++) {
    const Dual2 u(b0, 1.0, 0.0), w(b2, 0.0, 1.0);
    const Dual2 b[3] = {u, 1.0 - u - w, w};
    Dual2 X[3];
    kaseorg(b, X);
    const Dual2 f1 = X[i1] * c[m] - X[m] * c[i1];
    const Dual2 f2 = X[i2] * c[m] - X[m] * c[i2];
    const double det = f1.d[0] * f2.d[1] - f1.d[1] * f2.d[0];
    if (det == 0.0) break;
    const double du = (f1.v * f2.d[1] - f2.v * f1.d[1]) / det;
    const double dw = (f2.v * f1.d[0] - f1.v * f2.d[0]) / det;
    b0 -= du;
    b2 -= dw;
    if (std::fabs(du) + std::fabs(dw) < 1e-16) break;
  }
  return {0.5 * b0 + b2, kSin60 * b0};
}

// The lattice point under (z, az): AK's face point, through the warp's
// inverse when warped
std::pair<double,double> lattice_of_polar(double z, double az, bool warped) {
  const auto p = plane_of_polar(z, az);
  if (!warped) return p;
  double x, y;
  hex9::warp_solve(p.first, p.second, x, y);
  return {x, y};
}

} // anon

template <class T>
std::pair<T,T> ak_face_xy(T z, T az, bool warped) {
  if constexpr (std::is_same<T, double>::value) {
    return lattice_of_polar(z, az, warped);
  } else {
    // The derivatives of the forward are the inverse of those of the
    // plane-to-sphere direction: AK's in closed form, after the warp's
    // (I + grad d, read by central differences) when warped.
    const auto xy = lattice_of_polar(value_of(z), value_of(az), warped);
    double raw[2] = {xy.first, xy.second};
    double W[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
    if (warped) {
      constexpr double h = 1e-7;
      hex9::warp_apply(xy.first, xy.second, raw[0], raw[1]);
      double px, py, mx, my;
      hex9::warp_apply(xy.first + h, xy.second, px, py);
      hex9::warp_apply(xy.first - h, xy.second, mx, my);
      W[0][0] = (px - mx) / (2.0 * h);
      W[1][0] = (py - my) / (2.0 * h);
      hex9::warp_apply(xy.first, xy.second + h, px, py);
      hex9::warp_apply(xy.first, xy.second - h, mx, my);
      W[0][1] = (px - mx) / (2.0 * h);
      W[1][1] = (py - my) / (2.0 * h);
    }
    const auto pol = polar_of_plane(Dual2(raw[0], 1.0, 0.0), Dual2(raw[1], 0.0, 1.0));
    // d(z, az) / d(lattice) = d(z, az) / d(raw) . W
    const double a = pol.first.d[0] * W[0][0] + pol.first.d[1] * W[1][0];
    const double b = pol.first.d[0] * W[0][1] + pol.first.d[1] * W[1][1];
    const double c = pol.second.d[0] * W[0][0] + pol.second.d[1] * W[1][0];
    const double d = pol.second.d[0] * W[0][1] + pol.second.d[1] * W[1][1];
    const double det = a * d - b * c;
    T x, y;
    x.v = xy.first;
    y.v = xy.second;
    for (int k = 0; k < 2; k++) {
      x.d[k] = (d * z.d[k] - b * az.d[k]) / det;
      y.d[k] = (-c * z.d[k] + a * az.d[k]) / det;
    }
    return {x, y};
  }
}

template std::pair<double,double> ak_face_xy<double>(double, double, bool);
template std::pair<Dual2,Dual2> ak_face_xy<Dual2>(Dual2, Dual2, bool);

std::pair<double,double> ak_face_polar(double x, double y, bool warped) {
  if (warped) hex9::warp_apply(x, y, x, y);
  auto pol = polar_of_plane(x, y);
  if (pol.second < 0.0) pol.second += kTwoPi;
  return pol;
}

} // namespace hexify
