// hex_common.h - Shared utilities for hexagonal grid computations
//
// This file contains common constants, Class I quantization, and hex coordinate
// conversion functions used across aperture 3, 4, and 7 implementations.

#ifndef HEXIFY_HEX_COMMON_H
#define HEXIFY_HEX_COMMON_H

#include <cmath>

namespace hexify {
namespace common {

// =============================================================================
// Mathematical Constants
// =============================================================================

constexpr double M_SIN60 = 0.86602540378443864676;      // sin(60°)
constexpr double M_SQRT3 = 1.7320508075688772935;       // √3
constexpr double M_ONE_THIRD = 1.0 / 3.0;
constexpr double M_TWO_THIRDS = 2.0 / 3.0;
constexpr double PI = 3.14159265358979323846;

// =============================================================================
// Class I Hexagonal Grid Quantization (DGGRID DgHexC1Grid2D::quantify)
// =============================================================================

// Quantify a point (x, y) to Class I hexagonal grid cell (i, j)
// This is the core DGGRID algorithm for flat-top hex quantization
inline void quantify_class1_unit(double x, double y, long long& out_i, long long& out_j) {
  double a1 = std::fabs(x);
  double a2 = std::fabs(y);
  double x2 = a2 / M_SIN60;
  double x1 = a1 + x2 / 2.0;

  long long m1 = static_cast<long long>(x1);
  long long m2 = static_cast<long long>(x2);
  double r1 = x1 - m1;
  double r2 = x2 - m2;

  long long i_candidate, j_candidate;

  if (r1 < 0.5) {
    if (r1 < M_ONE_THIRD) {
      if (r2 < (1.0 + r1) / 2.0) {
        i_candidate = m1; j_candidate = m2;
      } else {
        i_candidate = m1; j_candidate = m2 + 1;
      }
    } else {
      j_candidate = (r2 < (1.0 - r1)) ? m2 : m2 + 1;
      i_candidate = ((1.0 - r1) <= r2 && r2 < (2.0 * r1)) ? m1 + 1 : m1;
    }
  } else {
    if (r1 < M_TWO_THIRDS) {
      j_candidate = (r2 < (1.0 - r1)) ? m2 : m2 + 1;
      i_candidate = ((2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)) ? m1 : m1 + 1;
    } else {
      if (r2 < (r1 / 2.0)) {
        i_candidate = m1 + 1; j_candidate = m2;
      } else {
        i_candidate = m1 + 1; j_candidate = m2 + 1;
      }
    }
  }

  // Handle negative coordinates
  if (x < 0.0) {
    if ((j_candidate % 2) == 0) {
      long long axisi = j_candidate / 2;
      i_candidate -= 2 * (i_candidate - axisi);
    } else {
      long long axisi = (j_candidate + 1) / 2;
      i_candidate -= 2 * (i_candidate - axisi) + 1;
    }
  }

  if (y < 0.0) {
    i_candidate -= (2 * j_candidate + 1) / 2;
    j_candidate = -j_candidate;
  }

  out_i = i_candidate;
  out_j = j_candidate;
}

// =============================================================================
// Class I Hex Center (DGGRID DgHexC1Grid2D::invQuantify)
// =============================================================================

// Convert hex cell (i, j) to center coordinates (x, y) in Class I grid
inline void ij_to_xy_class1(long long i, long long j, double& x, double& y) {
  x = static_cast<double>(i) - static_cast<double>(j) * 0.5;
  y = static_cast<double>(j) * M_SIN60;
}

// =============================================================================
// Class II Hexagonal Grid Quantization (DGGRID DgHexC2Grid2D::quantify)
// =============================================================================

// Quantify using Class II approach: rotate -30°, quantify in surrogate Class I,
// convert to substrate coordinates
inline void quantify_class2_unit(double x, double y, long long& out_i, long long& out_j) {
  // Step 1: Rotate by -30° (surrogate grid)
  const double angle = -30.0 * PI / 180.0;
  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);

  double x_rot = x * cos_a - y * sin_a;
  double y_rot = x * sin_a + y * cos_a;

  // Step 2: Quantize in surrogate Class I grid
  long long surr_i, surr_j;
  quantify_class1_unit(x_rot, y_rot, surr_i, surr_j);

  // Step 3: Convert surrogate to substrate (Class II output uses same i,j interpretation)
  // The substrate is at √3 finer scale, so surrogate coords map directly
  out_i = surr_i;
  out_j = surr_j;
}

// =============================================================================
// Hex Corners Generation
// =============================================================================

// Generate 6 corners of a hexagon centered at (cx, cy)
// rotation_offset: angle offset in radians (0 for Class I, PI/6 for Class II)
inline void generate_hex_corners(double cx, double cy, double hex_radius,
                                  double rotation_offset, double* xs, double* ys) {
  for (int k = 0; k < 6; ++k) {
    double angle = rotation_offset + k * PI / 3.0;
    xs[k] = cx + hex_radius * std::cos(angle);
    ys[k] = cy + hex_radius * std::sin(angle);
  }
}

} // namespace common
} // namespace hexify

#endif // HEXIFY_HEX_COMMON_H
