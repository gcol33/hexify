#include "authalic.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// Constants of the ellipsoid of flattening f.
struct Shape {
  double e2, e, qp;
  explicit Shape(double f) {
    if (!(f > 0.0 && f < 1.0)) {
      throw std::invalid_argument("ellipsoid flattening must lie in (0, 1)");
    }
    e2 = f * (2.0 - f);
    e = std::sqrt(e2);
    qp = 1.0 + (1.0 - e2) * std::atanh(e) / e;
  }
};

// atanh(x) / x, 1 at x = 0
inline double atanh_ratio(double x) {
  return x == 0.0 ? 1.0 : std::atanh(x) / x;
}

// q(phi) and r = (q_p - q(phi)) / cos^2(phi) at s = sin(phi) >= 0, c =
// cos(phi). Near the pole q_p - q is taken in closed form rather than as a
// difference, since atanh(e) - atanh(e s) = atanh(e (1 - s) / (1 - e^2 s)) and
// 1 - s = c^2 / (1 + s); r stays finite at the pole.
struct QTerms { double q, r; };

inline QTerms q_terms(const Shape& sh, double s, double c) {
  QTerms t;
  const double one_es2 = 1.0 - sh.e2 * s * s;
  if (s <= 0.5) {
    t.q = (1.0 - sh.e2) * (s / one_es2 + s * atanh_ratio(sh.e * s));
    t.r = (sh.qp - t.q) / (c * c);
  } else {
    const double denom = (1.0 + s) * (1.0 - sh.e2 * s);
    t.r = (1.0 + sh.e2 * s) / ((1.0 + s) * one_es2) +
          (1.0 - sh.e2) * atanh_ratio(sh.e * c * c / denom) / denom;
    t.q = sh.qp - t.r * c * c;
  }
  return t;
}

// Authalic latitude of phi >= 0, and its derivative.
inline double authalic_core(const Shape& sh, double phi, double* derivative) {
  const double s = std::sin(phi), c = std::cos(phi);
  const QTerms t = q_terms(sh, s, c);
  // cos(beta) = c sqrt(r (q_p + q)) / q_p, sin(beta) = q / q_p
  const double w = std::sqrt(t.r * (sh.qp + t.q));
  if (derivative) {
    const double one_es2 = 1.0 - sh.e2 * s * s;
    *derivative = 2.0 * (1.0 - sh.e2) / (one_es2 * one_es2 * w);
  }
  return std::atan2(t.q, c * w);
}

// The series a latitude conversion is evaluated by: beta - phi (or phi -
// beta) as a sine series in 2 phi (or 2 beta), its coefficients read off the
// closed form by a discrete sine transform at kNodes - 1 latitudes in
// (0, 90) degrees. The difference is analytic and periodic, so the
// coefficients fall off geometrically (as n^k, n the third flattening) and
// the transform is exact to the rounding of the closed form, about 1e-16
// radians at each latitude; terms below kTailBound, two units in the last
// place of one radian, are dropped.
constexpr int kNodes = 128;
constexpr int kMaxTerms = kNodes / 2;
constexpr double kTailBound = 4.4e-16;

struct LatitudeSeries {
  double flattening = 0.0;
  int n_forward = 0, n_inverse = 0;
  double forward[kMaxTerms];
  double inverse[kMaxTerms];
};

LatitudeSeries g_series;
bool g_active = false;

// Coefficients c_1..c_K of g(x) = sum_k c_k sin(2 k x) from g at x_j = j pi /
// (2 kNodes), j = 1..kNodes - 1; returns K.
int sine_coefficients(const double* g, double* out) {
  double sines[2 * kNodes];
  for (int m = 0; m < 2 * kNodes; ++m) sines[m] = std::sin(kPi * m / kNodes);
  int last = 0;
  for (int k = 1; k <= kMaxTerms; ++k) {
    double sum = 0.0;
    for (int j = 1; j < kNodes; ++j) sum += g[j] * sines[(k * j) % (2 * kNodes)];
    out[k - 1] = 2.0 * sum / kNodes;
    if (std::fabs(out[k - 1]) > kTailBound) last = k;
  }
  if (last == kMaxTerms) {
    throw std::invalid_argument(
        "ellipsoid flattening too large: its latitude series does not converge "
        "within 64 terms");
  }
  return last;
}

void build_series(double f) {
  double fwd[kNodes], inv[kNodes];
  for (int j = 1; j < kNodes; ++j) {
    const double x = j * kPiOver2 / kNodes;
    fwd[j] = authalic_lat_exact(x, f) - x;
    inv[j] = geodetic_lat_exact(x, f) - x;
  }
  LatitudeSeries s;
  s.n_forward = sine_coefficients(fwd, s.forward);
  s.n_inverse = sine_coefficients(inv, s.inverse);
  s.flattening = f;
  g_series = s;
}

// x + sum_k c_k sin(2 k x), by Clenshaw's recurrence in cos(2 x)
inline double sine_series(double x, const double* c, int n) {
  const double s = std::sin(x), co = std::cos(x);
  const double X = 2.0 * (co - s) * (co + s);
  double b1 = 0.0, b2 = 0.0;
  for (int k = n - 1; k >= 0; --k) {
    const double t = c[k] + X * b1 - b2;
    b2 = b1;
    b1 = t;
  }
  return x + 2.0 * s * co * b1;
}

} // namespace

void use_ellipsoid(double flattening) {
  if (flattening == 0.0) {
    g_active = false;
    return;
  }
  if (flattening != g_series.flattening) build_series(flattening);
  g_active = true;
}

double active_flattening() { return g_active ? g_series.flattening : 0.0; }

double to_sphere_lat(double lat_rad) {
  if (!g_active) return lat_rad;
  return sine_series(lat_rad, g_series.forward, g_series.n_forward);
}

double to_geodetic_lat(double lat_rad) {
  if (!g_active) return lat_rad;
  return sine_series(lat_rad, g_series.inverse, g_series.n_inverse);
}

double to_sphere_lat_deg(double lat_deg) {
  if (!g_active) return lat_deg;
  return to_sphere_lat(lat_deg * kDegToRad) * kRadToDeg;
}

double to_geodetic_lat_deg(double lat_deg) {
  if (!g_active) return lat_deg;
  return to_geodetic_lat(lat_deg * kDegToRad) * kRadToDeg;
}

double authalic_lat_exact(double phi, double flattening) {
  const Shape sh(flattening);
  const double beta = authalic_core(sh, std::fabs(phi), nullptr);
  return phi < 0.0 ? -beta : beta;
}

double geodetic_lat_exact(double beta, double flattening) {
  const Shape sh(flattening);
  const double b = std::fabs(beta);
  if (b >= kPiOver2) return beta;
  // Start from the series in e^2 to third order (PROJ's pj_authlat), then
  // Newton's method on the closed form.
  const double e2 = sh.e2, e4 = e2 * e2, e6 = e4 * e2;
  double phi = b + (e2 / 3.0 + 31.0 * e4 / 180.0 + 517.0 * e6 / 5040.0) * std::sin(2.0 * b) +
               (23.0 * e4 / 360.0 + 251.0 * e6 / 3780.0) * std::sin(4.0 * b) +
               (761.0 * e6 / 45360.0) * std::sin(6.0 * b);
  for (int iter = 0; iter < 60; ++iter) {
    double d;
    const double step = (authalic_core(sh, phi, &d) - b) / d;
    phi -= step;
    if (phi > kPiOver2) phi = kPiOver2;
    if (phi < 0.0) phi = 0.0;
    if (std::fabs(step) <= 4e-16) break;
  }
  return beta < 0.0 ? -phi : phi;
}

double authalic_radius_ratio(double flattening) {
  if (flattening == 0.0) return 1.0;
  return std::sqrt(Shape(flattening).qp / 2.0);
}

double authalic_parallel_scale(double phi, double flattening) {
  if (flattening == 0.0) return 1.0;
  const Shape sh(flattening);
  const double s = std::fabs(std::sin(phi)), c = std::cos(phi);
  const QTerms t = q_terms(sh, s, c);
  // k = R_q cos(beta) / (N cos(phi)), N = a / sqrt(1 - e^2 sin^2(phi))
  return std::sqrt(sh.qp / 2.0) * std::sqrt(1.0 - sh.e2 * s * s) *
         std::sqrt(t.r * (sh.qp + t.q)) / sh.qp;
}

} // namespace hexify
