// hex_ap34.cpp - Mixed aperture-3/4 hierarchical hexagon grid
//
// Handles arbitrary sequences of aperture 3 and 4 resolutions.
// Implements mixed aperture grids (e.g., ISEA43H) which use aperture 4 for coarse
// resolutions then switch to aperture 3 for finer resolutions.
//
// Hex orientations:
//   - Aperture 4: Always Rotation Class I (0-degree, flat-top)
//   - Aperture 3: Alternates Rotation Class I/II based on cumulative aperture-3 count
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "aperture_sequence.h"
#include "grid_math.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// ============================================================================
// Helper Functions
// ============================================================================

// Determine if a resolution uses 0-degree (flat-top) or 30-degree (pointy-top)
bool is_0deg_orientation(const std::vector<int>& ap_seq, size_t res_idx) {
    if (res_idx >= ap_seq.size()) {
        throw std::runtime_error("hex_ap34: resolution index out of bounds");
    }

    int current_ap = ap_seq[res_idx];

    if (current_ap == 4) {
        // Aperture 4 always uses 0-degree
        return true;
    } else if (current_ap == 3) {
        // Count aperture-3 resolutions up to and including this one
        int ap3_count = 0;
        for (size_t i = 0; i <= res_idx; i++) {
            if (ap_seq[i] == 3) {
                ap3_count++;
            }
        }
        // Odd count = 0-degree, even count = 30-degree
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

    // 30-degree orientation uses sqrt(3)x finer substrate
    size_t final_idx = ap_seq.size() - 1;
    if (!is_0deg_orientation(ap_seq, final_idx)) {
        scale *= kSqrt3;
    }

    return scale;
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

void hex_quantize_ap34(double icosa_triangle_x, double icosa_triangle_y,
                       const std::vector<int>& ap_seq,
                       long long& out_i, long long& out_j) {
    if (ap_seq.empty()) {
        throw std::runtime_error("hex_quantize_ap34: ap_seq cannot be empty");
    }

    double scale = calc_cumulative_scale(ap_seq);
    double grid_x = icosa_triangle_x * scale;
    double grid_y = icosa_triangle_y * scale;

    size_t final_idx = ap_seq.size() - 1;
    bool use_0deg = is_0deg_orientation(ap_seq, final_idx);

    if (use_0deg) {
        quantize_rotation_classI(grid_x, grid_y, out_i, out_j);
    } else {
        quantize_rotation_classII(grid_x, grid_y, out_i, out_j);
    }
}

void hex_center_ap34(long long i, long long j,
                     const std::vector<int>& ap_seq,
                     double& out_cx, double& out_cy) {
    if (ap_seq.empty()) {
        throw std::runtime_error("hex_center_ap34: ap_seq cannot be empty");
    }

    // (i,j) are in substrate coordinates (always Rotation Class I layout)
    center_rotation_classI(i, j, out_cx, out_cy);

    double substrate_scale = calc_substrate_scale(ap_seq);
    out_cx /= substrate_scale;
    out_cy /= substrate_scale;
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
    bool use_rotation_classI = is_0deg_orientation(ap_seq, final_idx);

    // Rotation Class I: 0 deg, Rotation Class II: 30 deg
    double rotation_deg = use_rotation_classI ? 0.0 : 30.0;

    generate_hex_corners(cx, cy, scaled_radius, rotation_deg, out_x, out_y);
}

} // namespace hexify
