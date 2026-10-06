#include "projection_fuller.h"
#include "constants.h"
#include <algorithm>
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

namespace {

// Gray's equation (39) in Crider's form (27): with u = a2 - alpha,
//   tan(u - A) + tan(u) + tan(u - B) + tan(alpha) = 0,
// A = x'' - sqrt(3) y'', B = 2 x''.
struct FullerShifts {
  double A, B;
};

FullerShifts fuller_shifts(double x, double y) {
  const double xpp = (x - ORIGIN_X) * kFullerArc;
  const double ypp = (y - ORIGIN_Y) * kFullerArc;
  return {xpp - SQRT3 * ypp, 2.0 * xpp};
}

// From t1 = tan(u - A) = tan(a1 - alpha) (Crider eq. 35) and t2 = tan(u) back
// to the secondary plane triangle (Crider eqs. 36-38), then onto the sphere
// along the ray through the plane point (Crider eq. 39).
std::pair<double,double> fuller_polar_from_tans(double t1, double t2) {
  const double a1s = kFullerDVE * t1;
  const double a2s = kFullerDVE * t2;
  const double ys = 0.5 * SQRT3 * (a1s + kFullerEL / 2.0 - kFullerEL / 3.0);
  const double xs = a2s + kFullerEL / 2.0 + ys / SQRT3 - kFullerEL / 3.0;

  const double z = std::atan2(std::hypot(xs, ys), kFullerZ0);
  double az = std::atan2(xs, ys);
  if (az < 0.0) az += kTwoPi;
  return {z, az};
}

// The closed-form inverse works in m = (b1 + b2 + b3) / 3, b_i = a_i - alpha,
// through t = tan(m - kFullerM0). Across the face, and on points up to a
// quarter of the face's size outside it, the face's root lies in
// m in [-0.204, -0.172] and the cubic's other two roots at |m| >= 1.09 (mod
// pi), so placing the point at infinity at m = 0.45 keeps the cubic's leading
// coefficient away from zero and makes the face's root its largest.
const double kFullerM0 = 0.45 - 0.5 * kPi;
const double kFullerT0 = std::tan(kFullerM0);
// (3 + ic) e^{3 i m0}, c = tan(alpha)
const double kFullerAr = 3.0 * std::cos(3.0 * kFullerM0) - kFullerTanAlpha * std::sin(3.0 * kFullerM0);
const double kFullerAi = 3.0 * std::sin(3.0 * kFullerM0) + kFullerTanAlpha * std::cos(3.0 * kFullerM0);
// (1 + ic) e^{i m0}
const double kFullerBr = std::cos(kFullerM0) - kFullerTanAlpha * std::sin(kFullerM0);
const double kFullerBi = std::sin(kFullerM0) + kFullerTanAlpha * std::cos(kFullerM0);

} // anon

std::pair<double,double> fuller_face_polar(double x, double y) {
  const auto [A, B] = fuller_shifts(x, y);

  // With b_i = m + d_i, d = (S/3 - A, S/3, S/3 - B), S = A + B, the equation
  // multiplied by 4 cos(b1) cos(b2) cos(b3) reads
  //   Im[(3 + ic) e^{3im} + (1 + ic) Z e^{im}] = 0,   Z = sum_i e^{-2i d_i}.
  // With tau_i = tan(d_i) and d3 = -(d1 + d2),
  //   e^{-2i d_1} = (1 - i tau1)^2 / s1,  e^{-2i d_2} = (1 - i tau2)^2 / s2,
  //   e^{-2i d_3} = (1 + i tau1)^2 (1 + i tau2)^2 / (s1 s2),  s_i = 1 + tau_i^2,
  // and the common factor s1 s2 is carried through every coefficient.
  const double S3 = (A + B) * (1.0 / 3.0);
  const double tau1 = std::tan(S3 - A);
  const double tau2 = std::tan(S3);
  const double s1 = 1.0 + tau1 * tau1;
  const double s2 = 1.0 + tau2 * tau2;
  const double w1r = 1.0 - tau1 * tau1, w1i = -2.0 * tau1;
  const double w2r = 1.0 - tau2 * tau2, w2i = -2.0 * tau2;
  const double w3r = w1r * w2r - w1i * w2i;
  const double w3i = -(w1r * w2i + w1i * w2r);
  const double zr = w1r * s2 + w2r * s1 + w3r;
  const double zi = w1i * s2 + w2i * s1 + w3i;
  const double den = s1 * s2;

  // In theta = m - m0, with a = (3 + ic) e^{3 i m0} and b = (1 + ic) Z e^{i m0},
  // the equation divided by cos^3(theta) is the cubic in t = tan(theta)
  //   (b_r - a_r) t^3 + (b_i - 3 a_i) t^2 + (3 a_r + b_r) t + (a_i + b_i) = 0.
  const double br = kFullerBr * zr - kFullerBi * zi;
  const double bi = kFullerBr * zi + kFullerBi * zr;
  const double ar = kFullerAr * den;
  const double ai = kFullerAi * den;
  const double c3 = br - ar;
  const double c2 = bi - 3.0 * ai;
  const double c1 = 3.0 * ar + br;
  const double c0 = ai + bi;

  // Its three roots are real; the largest, by Viete's trigonometric solution,
  //   t = (2 sgn(c3) sqrt(M) cos(phi / 3) - c2) / (3 c3),
  //   cos(phi) = -sgn(c3) N / (2 M^{3/2}),
  //   M = c2^2 - 3 c1 c3,  N = 2 c2^3 - 9 c1 c2 c3 + 27 c0 c3^2,
  // with cos(phi / 3) from tan(phi / 6). Where the other two roots meet, the
  // largest is stationary in phi, so it stays accurate there.
  const double M = c2 * c2 - 3.0 * c1 * c3;
  const double N = c2 * (2.0 * c2 * c2 - 9.0 * c1 * c3) + 27.0 * c0 * c3 * c3;
  const double sM = std::copysign(std::sqrt(M), c3);
  const double cos_phi = std::max(-1.0, std::min(1.0, -0.5 * N / (M * sM)));
  const double tt = std::tan(std::acos(cos_phi) * (1.0 / 6.0));
  const double tt2 = tt * tt;
  // t = tn / td
  const double tn = 2.0 * sM * (1.0 - tt2) - c2 * (1.0 + tt2);
  const double td = 3.0 * c3 * (1.0 + tt2);

  // tan(b_i) = tan((m0 + d_i) + theta), with tan(m0 + d_i) = n_i / e_i
  const double n1 = kFullerT0 + tau1, e1 = 1.0 - kFullerT0 * tau1;
  const double n2 = kFullerT0 + tau2, e2 = 1.0 - kFullerT0 * tau2;
  return fuller_polar_from_tans((n1 * td + tn * e1) / (e1 * td - tn * n1),
                                (n2 * td + tn * e2) / (e2 * td - tn * n2));
}

std::pair<double,double> fuller_face_polar_newton(double x, double y) {
  // The step after one of size d is of size ~d^2, so a step below 1e-14
  // leaves u at full double precision.
  constexpr double kTol = 1e-14;
  constexpr int kMaxIters = 100;

  // Every term increases with u, so the root is unique and Newton's method
  // converges from the root of the linearised equation.
  const auto [A, B] = fuller_shifts(x, y);
  double u = (A + B - kFullerTanAlpha) / 3.0;
  for (int k = 0; k < kMaxIters; ++k) {
    const double t1 = std::tan(u - A);
    const double t2 = std::tan(u);
    const double t3 = std::tan(u - B);
    const double g = t1 + t2 + t3 + kFullerTanAlpha;
    const double dg = 3.0 + t1 * t1 + t2 * t2 + t3 * t3;
    const double du = g / dg;
    u -= du;
    if (std::fabs(du) <= kTol) break;
  }
  return fuller_polar_from_tans(std::tan(u - A), std::tan(u));
}

} // namespace hexify
