// hex_ap7.cpp - Aperture-7 hierarchical hexagon grid
// Based on DGGRID DgHexC3Grid2D: uses Class III with alternating variants
// Class III-I (even res): substrate is sqrt(7)× finer
// Class III-II (odd res): substrate is sqrt(7)×sqrt(3)× finer = sqrt(21)× finer

#include "hex_ap7.h"
#include <cmath>
#include <stdexcept>

namespace {
  const double M_SIN60 = 0.86602540378443864676;
  const double M_SQRT3 = 1.7320508075688772935;
  const double M_SQRT7 = 2.6457513110645905905;
  const double M_SQRT21 = 4.5825756949558400065;  // sqrt(7*3)
  const double M_ONE_THIRD = 1.0 / 3.0;
  const double M_TWO_THIRDS = 2.0 / 3.0;
  const double PI = 3.14159265358979323846;
  
  // M_AP7_ROT_DEGS from DGGRID - need to find this value
  // From geometry: arctan(sqrt(3/7)) ≈ 19.1066° 
  const double M_AP7_ROT_DEG = 19.10660535;
}

namespace hexify {

// ============== Class I Quantify (used as substrate) ==============

static void quantify_class1_unit(double x, double y, long long& out_i, long long& out_j) {
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

static void ij_to_xy_class1(long long i, long long j, double& x, double& y) {
  x = static_cast<double>(i) - static_cast<double>(j) * 0.5;
  y = static_cast<double>(j) * M_SIN60;
}

// ============== Class II Quantify (used for Class III-II surrogate) ==============

static void quantify_class2_unit(double x, double y, long long& out_i, long long& out_j) {
  const double angle = -30.0 * PI / 180.0;
  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);
  
  double rot_x = x * cos_a - y * sin_a;
  double rot_y = x * sin_a + y * cos_a;
  
  long long sur_i, sur_j;
  quantify_class1_unit(rot_x, rot_y, sur_i, sur_j);
  
  double sur_cx, sur_cy;
  ij_to_xy_class1(sur_i, sur_j, sur_cx, sur_cy);
  
  const double cos_back = std::cos(-angle);
  const double sin_back = std::sin(-angle);
  
  double back_x = sur_cx * cos_back - sur_cy * sin_back;
  double back_y = sur_cx * sin_back + sur_cy * cos_back;
  
  double sub_x = back_x * M_SQRT3;
  double sub_y = back_y * M_SQRT3;
  
  quantify_class1_unit(sub_x, sub_y, out_i, out_j);
}

// ============== Class III-I (Even Resolutions) ==============
// Surrogate: Class I rotated by M_AP7_ROT_DEG
// Substrate: Class I at sqrt(7)× finer scale

static void quantify_class3i_unit(double x, double y, long long& out_i, long long& out_j) {
  const double angle = -M_AP7_ROT_DEG * PI / 180.0;
  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);
  
  // Rotate to surrogate
  double rot_x = x * cos_a - y * sin_a;
  double rot_y = x * sin_a + y * cos_a;
  
  // Quantify in Class I surrogate
  long long sur_i, sur_j;
  quantify_class1_unit(rot_x, rot_y, sur_i, sur_j);
  
  // Get surrogate center
  double sur_cx, sur_cy;
  ij_to_xy_class1(sur_i, sur_j, sur_cx, sur_cy);
  
  // Rotate back
  const double cos_back = std::cos(-angle);
  const double sin_back = std::sin(-angle);
  
  double back_x = sur_cx * cos_back - sur_cy * sin_back;
  double back_y = sur_cx * sin_back + sur_cy * cos_back;
  
  // Scale to substrate (sqrt(7)× finer)
  double sub_x = back_x * M_SQRT7;
  double sub_y = back_y * M_SQRT7;
  
  // Quantify in substrate (Class I)
  quantify_class1_unit(sub_x, sub_y, out_i, out_j);
}

// ============== Class III-II (Odd Resolutions) ==============
// Surrogate: Class II rotated by M_AP7_ROT_DEG
// Substrate: Class I at sqrt(21)× finer scale = sqrt(7)×sqrt(3)×

static void quantify_class3ii_unit(double x, double y, long long& out_i, long long& out_j) {
  const double angle = -M_AP7_ROT_DEG * PI / 180.0;
  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);
  
  // Rotate to surrogate
  double rot_x = x * cos_a - y * sin_a;
  double rot_y = x * sin_a + y * cos_a;
  
  // Quantify in Class II surrogate
  long long sur_i, sur_j;
  quantify_class2_unit(rot_x, rot_y, sur_i, sur_j);
  
  // Get surrogate center (Class II, so substrate coords)
  double sur_cx, sur_cy;
  ij_to_xy_class1(sur_i, sur_j, sur_cx, sur_cy);
  
  // Surrogate center is in Class II substrate coords (sqrt(3)× finer than base)
  // Need to scale to our substrate coordinates
  // Our substrate should be sqrt(7)× finer than surrogate base
  // Since surrogate substrate is already sqrt(3)×, multiply by sqrt(7)
  
  // Rotate back
  const double cos_back = std::cos(-angle);
  const double sin_back = std::sin(-angle);
  
  double back_x = sur_cx * cos_back - sur_cy * sin_back;
  double back_y = sur_cx * sin_back + sur_cy * cos_back;
  
  // Scale to substrate (sqrt(7)× finer than Class II base)
  // Total: sqrt(3) * sqrt(7) = sqrt(21)
  double sub_x = back_x * M_SQRT7;
  double sub_y = back_y * M_SQRT7;
  
  // Quantify in substrate (Class I)
  quantify_class1_unit(sub_x, sub_y, out_i, out_j);
}

// ============== Public API ==============

void hex_quantify_ap7(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j) {
  if (resolution < 0) {
    throw std::runtime_error("hex_quantify_ap7: resolution must be >= 0");
  }
  
  // Cumulative scale: sqrt(7)^resolution
  double scale = std::pow(M_SQRT7, resolution);
  double scaled_x = tx * scale;
  double scaled_y = ty * scale;
  
  // Determine which Class III variant to use
  bool is_class3i = (resolution % 2 == 0);
  
  if (is_class3i) {
    quantify_class3i_unit(scaled_x, scaled_y, out_i, out_j);
  } else {
    quantify_class3ii_unit(scaled_x, scaled_y, out_i, out_j);
  }
}

void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
  if (resolution < 0) {
    throw std::runtime_error("hex_center_ap7: resolution must be >= 0");
  }
  
  // (i,j) are always in substrate coordinates (Class I)
  double x, y;
  ij_to_xy_class1(i, j, x, y);
  
  // The substrate scale is: base_scale × substrate_multiplier
  // - base_scale = sqrt(7)^resolution (from DgHexGrid2DS)
  // - substrate_multiplier depends on Class III variant (from DgHexC3Grid2D)
  
  double base_scale = std::pow(M_SQRT7, resolution);
  
  bool is_class3i = (resolution % 2 == 0);
  double substrate_multiplier;
  
  if (is_class3i) {
    // Class III-I: substrate is sqrt(7)× finer than base
    substrate_multiplier = M_SQRT7;
  } else {
    // Class III-II: substrate is sqrt(21)× finer than base
    substrate_multiplier = M_SQRT21;
  }
  
  double substrate_scale = base_scale * substrate_multiplier;
  
  out_cx = x / substrate_scale;
  out_cy = y / substrate_scale;
}

void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
  if (resolution < 0) {
    throw std::runtime_error("hex_corners_ap7: resolution must be >= 0");
  }
  
  double cx, cy;
  hex_center_ap7(i, j, resolution, cx, cy);
  
  double scale = std::pow(M_SQRT7, resolution);
  double scaled_radius = hex_radius / scale;
  
  // Class III hexagons are rotated by M_AP7_ROT_DEG
  // But we need to check which variant for proper orientation
  bool is_class3i = (resolution % 2 == 0);
  double rotation_offset = is_class3i ? M_AP7_ROT_DEG : (M_AP7_ROT_DEG + 30.0);
  
  for (int k = 0; k < 6; ++k) {
    double angle = (90.0 + rotation_offset + k * 60.0) * PI / 180.0;
    out_x[k] = cx + scaled_radius * std::cos(angle);
    out_y[k] = cy + scaled_radius * std::sin(angle);
  }
}

} // namespace hexify
