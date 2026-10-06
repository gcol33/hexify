#pragma once
// A number carrying its derivatives along two directions (forward-mode
// automatic differentiation). Code written for a scalar type T, run with
// T = Dual2, returns its value and its exact first derivatives: every
// operation applies its own derivative rule, so the result is the analytic
// Jacobian of the code as written, without a step size.
#include <cmath>

namespace hexify {

struct Dual2 {
  double v;      // value
  double d[2];   // derivatives along the two directions

  explicit Dual2(double x = 0.0) : v(x), d{0.0, 0.0} {}
  Dual2(double x, double d0, double d1) : v(x), d{d0, d1} {}
};

inline double value_of(double x) { return x; }
inline double value_of(const Dual2& x) { return x.v; }

// f(x) with f'(x) = df
inline Dual2 chain(const Dual2& x, double f, double df) {
  return {f, df * x.d[0], df * x.d[1]};
}

inline Dual2 operator+(const Dual2& a, const Dual2& b) {
  return {a.v + b.v, a.d[0] + b.d[0], a.d[1] + b.d[1]};
}
inline Dual2 operator-(const Dual2& a, const Dual2& b) {
  return {a.v - b.v, a.d[0] - b.d[0], a.d[1] - b.d[1]};
}
inline Dual2 operator-(const Dual2& a) { return {-a.v, -a.d[0], -a.d[1]}; }
inline Dual2 operator*(const Dual2& a, const Dual2& b) {
  return {a.v * b.v, a.d[0] * b.v + a.v * b.d[0], a.d[1] * b.v + a.v * b.d[1]};
}
inline Dual2 operator/(const Dual2& a, const Dual2& b) {
  const double q = a.v / b.v;
  return {q, (a.d[0] - q * b.d[0]) / b.v, (a.d[1] - q * b.d[1]) / b.v};
}
inline Dual2 operator+(const Dual2& a, double b) { return {a.v + b, a.d[0], a.d[1]}; }
inline Dual2 operator+(double a, const Dual2& b) { return b + a; }
inline Dual2 operator-(const Dual2& a, double b) { return {a.v - b, a.d[0], a.d[1]}; }
inline Dual2 operator-(double a, const Dual2& b) { return {a - b.v, -b.d[0], -b.d[1]}; }
inline Dual2 operator*(const Dual2& a, double b) { return {a.v * b, a.d[0] * b, a.d[1] * b}; }
inline Dual2 operator*(double a, const Dual2& b) { return b * a; }
inline Dual2 operator/(const Dual2& a, double b) { return {a.v / b, a.d[0] / b, a.d[1] / b}; }
inline Dual2 operator/(double a, const Dual2& b) { return Dual2(a) / b; }
inline Dual2& operator+=(Dual2& a, const Dual2& b) { return a = a + b; }
inline Dual2& operator+=(Dual2& a, double b) { return a = a + b; }

inline Dual2 sin(const Dual2& x) { return chain(x, std::sin(x.v), std::cos(x.v)); }
inline Dual2 cos(const Dual2& x) { return chain(x, std::cos(x.v), -std::sin(x.v)); }
inline Dual2 tan(const Dual2& x) {
  const double t = std::tan(x.v);
  return chain(x, t, 1.0 + t * t);
}
inline Dual2 atan(const Dual2& x) { return chain(x, std::atan(x.v), 1.0 / (1.0 + x.v * x.v)); }
inline Dual2 acos(const Dual2& x) {
  return chain(x, std::acos(x.v), -1.0 / std::sqrt(1.0 - x.v * x.v));
}
inline Dual2 sqrt(const Dual2& x) {
  const double s = std::sqrt(x.v);
  return chain(x, s, 0.5 / s);
}
inline Dual2 atan2(const Dual2& y, const Dual2& x) {
  const double r2 = x.v * x.v + y.v * y.v;
  return {std::atan2(y.v, x.v), (x.v * y.d[0] - y.v * x.d[0]) / r2,
          (x.v * y.d[1] - y.v * x.d[1]) / r2};
}
inline Dual2 atan2(double y, const Dual2& x) { return atan2(Dual2(y), x); }

// x held to [lo, hi]; held at a bound, it no longer moves.
inline double clamp_to(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}
inline Dual2 clamp_to(const Dual2& x, double lo, double hi) {
  if (x.v < lo) return Dual2(lo);
  if (x.v > hi) return Dual2(hi);
  return x;
}

} // namespace hexify
