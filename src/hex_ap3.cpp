// hex_ap3.cpp - Aperture-3 hierarchical hexagon grid
// Pure DGGRID port: works with (i,j) coordinates at each resolution
// Includes both Class I (even resolutions) and Class II (odd resolutions)
#include "hex_ap3.h"
#include <cmath>
#include <stdexcept>

namespace {
  const double M_SIN60 = 0.86602540378443864676;
  const double M_SQRT3 = 1.7320508075688772935;
  const double M_ONE_THIRD = 1.0 / 3.0;
  const double M_TWO_THIRDS = 2.0 / 3.0;
  const double PI = 3.14159265358979323846;
}

namespace hexify {

// ============== Class I (Even Resolutions) ==============

// DGGRID's DgHexC1Grid2D::quantify() - exact port
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

// DGGRID's coordinate conversion for Class I hex
// Center of hex (i,j) in unit grid
static void ij_to_xy_class1(long long i, long long j, double& x, double& y) {
  x = static_cast<double>(i) - static_cast<double>(j) * 0.5;  // FIXED: subtract, not add
  y = static_cast<double>(j) * M_SIN60;
}

// ============== Class II (Odd Resolutions) ==============

// DGGRID's DgHexC2Grid2D approach:
// Class II uses a "surrogate" Class I grid rotated -30 degrees for quantization
// Then converts to "substrate" Class I grid (one aperture 3 finer = scaled by √3)
//
// Key insight from DGGRID source:
// - Surrogate: Class I at same scale, rotated -30°
// - Substrate: Class I at √3 finer scale (one resolution higher)
// - Quantify: rotate → quantify in surrogate → express in substrate coords

static void quantify_class2_unit(double x, double y, long long& out_i, long long& out_j) {
  // DGGRID Class II quantify: rotate -30° → quantify in Class I → convert to substrate coords
  
  // Step 1: Rotate by -30° (surrogate grid)
  const double angle = -30.0 * PI / 180.0;
  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);
  
  double rot_x = x * cos_a - y * sin_a;
  double rot_y = x * sin_a + y * cos_a;
  
  // Step 2: Quantify in the surrogate (Class I rotated) grid
  long long sur_i, sur_j;
  quantify_class1_unit(rot_x, rot_y, sur_i, sur_j);
  
  // Step 3: Get surrogate center in world coords
  double sur_cx, sur_cy;
  ij_to_xy_class1(sur_i, sur_j, sur_cx, sur_cy);
  
  // Step 4: Rotate surrogate center back (+30°)
  const double cos_back = std::cos(-angle);
  const double sin_back = std::sin(-angle);
  
  double back_x = sur_cx * cos_back - sur_cy * sin_back;
  double back_y = sur_cx * sin_back + sur_cy * cos_back;
  
  // Step 5: Express in substrate coordinates (Class I at √3 finer scale)
  double sub_x = back_x * M_SQRT3;
  double sub_y = back_y * M_SQRT3;
  
  // Step 6: Quantify in substrate to get final address
  quantify_class1_unit(sub_x, sub_y, out_i, out_j);
}

static void ij_to_xy_class2(long long i, long long j, double& x, double& y) {
  // DGGRID Class II invQuantify: (i,j) are substrate coordinates
  // Simply get Class I position at substrate scale (which is √3 finer)
  // The substrate IS a Class I grid, so just use Class I formula
  ij_to_xy_class1(i, j, x, y);
}

// ============== Public API ==============

void hex_quantify_ap3(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j) {
  if (resolution < 0) {
    throw std::runtime_error("hex_quantify_ap3: resolution must be >= 0");
  }
  
  // Determine class based on resolution
  bool is_class1 = (resolution % 2 == 0);
  
  // DGGRID scaling: each resolution r has cumulative scale sqrt(3)^r
  double scale = std::pow(M_SQRT3, resolution);
  
  // Scale point to resolution grid
  double scaled_x = tx * scale;
  double scaled_y = ty * scale;
  
  // Quantify using appropriate class
  if (is_class1) {
    // Even resolution: Class I
    quantify_class1_unit(scaled_x, scaled_y, out_i, out_j);
  } else {
    // Odd resolution: Class II
    quantify_class2_unit(scaled_x, scaled_y, out_i, out_j);
  }
}

void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
  if (resolution < 0) {
    throw std::runtime_error("hex_center_ap3: resolution must be >= 0");
  }
  
  // Determine class based on resolution
  bool is_class1 = (resolution % 2 == 0);
  
  // Get center in unit grid using appropriate class
  double x, y;
  if (is_class1) {
    // Class I: straightforward
    ij_to_xy_class1(i, j, x, y);
    
    // Scale back to base coordinates
    double scale = std::pow(M_SQRT3, resolution);
    out_cx = x / scale;
    out_cy = y / scale;
  } else {
    // Class II: (i,j) are in substrate coordinates (Class I at √3 finer scale)
    ij_to_xy_class2(i, j, x, y);
    
    // Substrate is √3 finer, so effective scale is sqrt(3)^(resolution+1)
    // But we can also think of it as: sqrt(3)^resolution * sqrt(3)
    double scale = std::pow(M_SQRT3, resolution) * M_SQRT3;
    out_cx = x / scale;
    out_cy = y / scale;
  }
}

void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
  if (resolution < 0) {
    throw std::runtime_error("hex_corners_ap3: resolution must be >= 0");
  }
  
  // Get center
  double cx, cy;
  hex_center_ap3(i, j, resolution, cx, cy);
  
  // DGGRID hex vertices
  // The hex_radius parameter is in the SCALED space (unit grid)
  // We need to scale it down to the resolution space
  
  double scale = std::pow(M_SQRT3, resolution);
  double scaled_radius = hex_radius / scale;
  
  // Class II hexagons are rotated 30° relative to Class I
  bool is_class1 = (resolution % 2 == 0);
  double rotation_offset = is_class1 ? 0.0 : 30.0;
  
  // Generate 6 vertices (flat-top orientation for Class I, rotated for Class II)
  // Starting at top, going counter-clockwise
  for (int k = 0; k < 6; ++k) {
    double angle = (90.0 + rotation_offset + k * 60.0) * PI / 180.0;
    out_x[k] = cx + scaled_radius * std::cos(angle);
    out_y[k] = cy + scaled_radius * std::sin(angle);
  }
}

} // namespace hexify
