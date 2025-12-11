// hex_ap34.cpp - Mixed aperture-3/4 hierarchical hexagon grid
// Handles arbitrary sequences of aperture 3 and 4 resolutions
// Following DGGRID's approach: ap4 uses Class I, ap3 alternates Class I/II

#include "hex_ap34.h"
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

// ============== Class I Quantify ==============

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

// ============== Class II Quantify ==============

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

// ============== Helper Functions ==============

// Determine if a resolution uses Class I or Class II
// Based on DGGRID's logic in DgHexGrid2DS.cpp
static bool is_class_one(const std::vector<int>& ap_seq, size_t res_idx) {
  if (res_idx >= ap_seq.size()) {
    throw std::runtime_error("hex_ap34: resolution index out of bounds");
  }
  
  int current_ap = ap_seq[res_idx];
  
  if (current_ap == 4) {
    // Aperture 4 always uses Class I
    return true;
  } else if (current_ap == 3) {
    // For aperture 3: count how many aperture-3 resolutions INCLUDING this one
    int ap3_count = 0;
    for (size_t i = 0; i <= res_idx; i++) {
      if (ap_seq[i] == 3) {
        ap3_count++;
      }
    }
    // First ap3 is Class I (count=1, odd), second is Class II (count=2, even), etc.
    // So: odd count = Class I, even count = Class II
    return (ap3_count % 2) == 1;
  } else {
    throw std::runtime_error("hex_ap34: aperture must be 3 or 4");
  }
}

// Calculate cumulative scale from aperture sequence
// The sequence represents apertures at each resolution 0, 1, 2, ...
// But resolution 0 is the base grid with no refinement
// So we only multiply by apertures starting from index 1
static double calc_cumulative_scale(const std::vector<int>& ap_seq) {
  double scale = 1.0;
  // Start from index 1, not 0! Resolution 0 is base grid.
  for (size_t i = 1; i < ap_seq.size(); i++) {
    int ap = ap_seq[i];
    if (ap == 3) {
      scale *= M_SQRT3;
    } else if (ap == 4) {
      scale *= 2.0;
    } else {
      throw std::runtime_error("hex_ap34: aperture must be 3 or 4");
    }
  }
  return scale;
}

// Calculate substrate scale at final resolution
static double calc_substrate_scale(const std::vector<int>& ap_seq) {
  if (ap_seq.empty()) {
    return 1.0;
  }
  
  double scale = calc_cumulative_scale(ap_seq);
  
  // If final resolution is Class II, substrate is sqrt(3)× finer
  size_t final_idx = ap_seq.size() - 1;
  if (!is_class_one(ap_seq, final_idx)) {
    scale *= M_SQRT3;
  }
  
  return scale;
}

// ============== Public API ==============

void hex_quantify_ap34(double tx, double ty, 
                       const std::vector<int>& ap_seq,
                       long long& out_i, long long& out_j) {
  if (ap_seq.empty()) {
    throw std::runtime_error("hex_quantify_ap34: ap_seq cannot be empty");
  }
  
  // Calculate cumulative scale
  double scale = calc_cumulative_scale(ap_seq);
  
  // Scale point
  double scaled_x = tx * scale;
  double scaled_y = ty * scale;
  
  // Determine which class to use at final resolution
  size_t final_idx = ap_seq.size() - 1;
  bool use_class1 = is_class_one(ap_seq, final_idx);
  
  // Quantify
  if (use_class1) {
    quantify_class1_unit(scaled_x, scaled_y, out_i, out_j);
  } else {
    quantify_class2_unit(scaled_x, scaled_y, out_i, out_j);
  }
}

void hex_center_ap34(long long i, long long j,
                     const std::vector<int>& ap_seq,
                     double& out_cx, double& out_cy) {
  if (ap_seq.empty()) {
    throw std::runtime_error("hex_center_ap34: ap_seq cannot be empty");
  }
  
  // (i,j) are in substrate coordinates (always Class I)
  double x, y;
  ij_to_xy_class1(i, j, x, y);
  
  // Scale back using substrate scale
  double substrate_scale = calc_substrate_scale(ap_seq);
  out_cx = x / substrate_scale;
  out_cy = y / substrate_scale;
}

void hex_corners_ap34(long long i, long long j,
                      const std::vector<int>& ap_seq,
                      double hex_radius,
                      double* out_x, double* out_y) {
  if (ap_seq.empty()) {
    throw std::runtime_error("hex_corners_ap34: ap_seq cannot be empty");
  }
  
  // Get center
  double cx, cy;
  hex_center_ap34(i, j, ap_seq, cx, cy);
  
  // Calculate scale for radius adjustment
  double scale = calc_cumulative_scale(ap_seq);
  double scaled_radius = hex_radius / scale;
  
  // Determine rotation based on final resolution class
  size_t final_idx = ap_seq.size() - 1;
  bool use_class1 = is_class_one(ap_seq, final_idx);
  double rotation_offset = use_class1 ? 0.0 : 30.0;
  
  // Generate 6 vertices
  for (int k = 0; k < 6; ++k) {
    double angle = (90.0 + rotation_offset + k * 60.0) * PI / 180.0;
    out_x[k] = cx + scaled_radius * std::cos(angle);
    out_y[k] = cy + scaled_radius * std::sin(angle);
  }
}

} // namespace hexify
