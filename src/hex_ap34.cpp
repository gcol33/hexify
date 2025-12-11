// hex_ap34.cpp - Mixed aperture-3/4 hierarchical hexagon grid
//
// Handles arbitrary sequences of aperture 3 and 4 resolutions.
// Implementation based on cube coordinate systems for hexagonal grids.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "hex_ap34.h"
#include "hex_coord.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// Class I (flat-top) quantization
void quantify_class1(double x, double y, long long& out_i, long long& out_j) {
    CubeCoord cube = cartesian_to_cube(x, y, kSqrt3);
    cube.round_to_nearest();
    out_i = static_cast<long long>(cube.q);
    out_j = static_cast<long long>(cube.r);
}

void center_class1(long long i, long long j, double& x, double& y) {
    cube_to_cartesian(static_cast<double>(i), static_cast<double>(j), x, y, kSin60);
}

// Class II (30° rotated) quantization
void quantify_class2(double x, double y, long long& out_i, long long& out_j) {
    constexpr double angle = -kPi / 6.0;  // -30°
    double c = std::cos(angle);
    double s = std::sin(angle);

    // Rotate to surrogate Class I
    double rx = x * c - y * s;
    double ry = x * s + y * c;

    // Quantize in surrogate
    long long sur_i, sur_j;
    quantify_class1(rx, ry, sur_i, sur_j);

    // Get surrogate center and rotate back
    double sur_x, sur_y;
    center_class1(sur_i, sur_j, sur_x, sur_y);

    double back_x = sur_x * c + sur_y * s;
    double back_y = -sur_x * s + sur_y * c;

    // Scale to substrate and re-quantize
    quantify_class1(back_x * kSqrt3, back_y * kSqrt3, out_i, out_j);
}

// ============================================================================
// Helper Functions
// ============================================================================

// Determine if a resolution uses Class I or Class II
bool is_class_one(const std::vector<int>& ap_seq, size_t res_idx) {
    if (res_idx >= ap_seq.size()) {
        throw std::runtime_error("hex_ap34: resolution index out of bounds");
    }

    int current_ap = ap_seq[res_idx];

    if (current_ap == 4) {
        // Aperture 4 always uses Class I
        return true;
    } else if (current_ap == 3) {
        // Count aperture-3 resolutions up to and including this one
        int ap3_count = 0;
        for (size_t i = 0; i <= res_idx; i++) {
            if (ap_seq[i] == 3) {
                ap3_count++;
            }
        }
        // Odd count = Class I, even count = Class II
        return (ap3_count % 2) == 1;
    } else {
        throw std::runtime_error("hex_ap34: aperture must be 3 or 4");
    }
}

// Calculate cumulative scale from aperture sequence
double calc_cumulative_scale(const std::vector<int>& ap_seq) {
    double scale = 1.0;
    // Start from index 1 (resolution 0 is base grid with no refinement)
    for (size_t i = 1; i < ap_seq.size(); i++) {
        int ap = ap_seq[i];
        if (ap == 3) {
            scale *= kSqrt3;
        } else if (ap == 4) {
            scale *= 2.0;
        } else {
            throw std::runtime_error("hex_ap34: aperture must be 3 or 4");
        }
    }
    return scale;
}

// Calculate substrate scale at final resolution
double calc_substrate_scale(const std::vector<int>& ap_seq) {
    if (ap_seq.empty()) {
        return 1.0;
    }

    double scale = calc_cumulative_scale(ap_seq);

    // Class II uses sqrt(3)× finer substrate
    size_t final_idx = ap_seq.size() - 1;
    if (!is_class_one(ap_seq, final_idx)) {
        scale *= kSqrt3;
    }

    return scale;
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

void hex_quantify_ap34(double tx, double ty,
                       const std::vector<int>& ap_seq,
                       long long& out_i, long long& out_j) {
    if (ap_seq.empty()) {
        throw std::runtime_error("hex_quantify_ap34: ap_seq cannot be empty");
    }

    double scale = calc_cumulative_scale(ap_seq);
    double scaled_x = tx * scale;
    double scaled_y = ty * scale;

    size_t final_idx = ap_seq.size() - 1;
    bool use_class1 = is_class_one(ap_seq, final_idx);

    if (use_class1) {
        quantify_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        quantify_class2(scaled_x, scaled_y, out_i, out_j);
    }
}

void hex_center_ap34(long long i, long long j,
                     const std::vector<int>& ap_seq,
                     double& out_cx, double& out_cy) {
    if (ap_seq.empty()) {
        throw std::runtime_error("hex_center_ap34: ap_seq cannot be empty");
    }

    // (i,j) are in substrate coordinates (always Class I layout)
    double x, y;
    center_class1(i, j, x, y);

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

    double cx, cy;
    hex_center_ap34(i, j, ap_seq, cx, cy);

    double scale = calc_cumulative_scale(ap_seq);
    double scaled_radius = hex_radius / scale;

    size_t final_idx = ap_seq.size() - 1;
    bool use_class1 = is_class_one(ap_seq, final_idx);
    double rotation_offset = use_class1 ? 0.0 : (kPi / 6.0);

    for (int k = 0; k < 6; ++k) {
        double angle = (kPi / 2.0) + rotation_offset + k * (kPi / 3.0);
        out_x[k] = cx + scaled_radius * std::cos(angle);
        out_y[k] = cy + scaled_radius * std::sin(angle);
    }
}

} // namespace hexify
