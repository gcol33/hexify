// hex_ap4.cpp - Aperture-4 hexagon grid
// Based on DGGRID DgHexGrid2DS.cpp - aperture 4 always uses Class I with scale factor 2

#include "hex_ap4.h"
#include <cmath>
#include <stdexcept>

namespace {
  const double M_SIN60 = 0.86602540378443864676;
  const double M_ONE_THIRD = 1.0 / 3.0;
  const double M_TWO_THIRDS = 2.0 / 3.0;
  const double PI = 3.14159265358979323846;
  
  // Aperture 4 scale factor (from DGGRID DgHexGrid2DS.cpp line 88)
  const double AP4_SCALE_FACTOR = 2.0;
}

namespace hexify {

// ============== Aperture 4 Grid ==============
// DGGRID insight: Aperture 4 ALWAYS uses Class I (no Class II)
// Scale factor: 2^resolution (not sqrt(3)^resolution)

// DGGRID's DgHexC1Grid2D::quantify() - same as aperture 3 Class I
static void quantify_ap4_unit(double x, double y, long long& out_i, long long& out_j) {
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

// DGGRID's DgHexC1Grid2D::invQuantify() - Class I formula
static void ij_to_xy_ap4(long long i, long long j, double& x, double& y) {
  x = static_cast<double>(i) - static_cast<double>(j) * 0.5;
  y = static_cast<double>(j) * M_SIN60;
}

// ============== Public API ==============

void hex_quantify_ap4(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j) {
  if (resolution < 0) {
    throw std::runtime_error("hex_quantify_ap4: resolution must be >= 0");
  }
  
  // Aperture 4: scale factor is 2^resolution (DGGRID DgHexGrid2DS.cpp line 88)
  double scale = std::pow(AP4_SCALE_FACTOR, resolution);
  
  double scaled_x = tx * scale;
  double scaled_y = ty * scale;
  
  // Quantify using Class I logic (aperture 4 is always Class I)
  quantify_ap4_unit(scaled_x, scaled_y, out_i, out_j);
}

void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
  if (resolution < 0) {
    throw std::runtime_error("hex_center_ap4: resolution must be >= 0");
  }
  
  // Get center in unit grid (Class I)
  double x, y;
  ij_to_xy_ap4(i, j, x, y);
  
  // Scale back: aperture 4 uses 2^resolution
  double scale = std::pow(AP4_SCALE_FACTOR, resolution);
  out_cx = x / scale;
  out_cy = y / scale;
}

void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
  if (resolution < 0) {
    throw std::runtime_error("hex_corners_ap4: resolution must be >= 0");
  }
  
  // Get center
  double cx, cy;
  hex_center_ap4(i, j, resolution, cx, cy);
  
  // Scale radius
  double scale = std::pow(AP4_SCALE_FACTOR, resolution);
  double scaled_radius = hex_radius / scale;
  
  // Flat-top hexagon (Class I orientation, same as aperture 3 Class I)
  // Starting at top, going counter-clockwise
  for (int k = 0; k < 6; ++k) {
    double angle = (90.0 + k * 60.0) * PI / 180.0;
    out_x[k] = cx + scaled_radius * std::cos(angle);
    out_y[k] = cy + scaled_radius * std::sin(angle);
  }
}

} // namespace hexify
