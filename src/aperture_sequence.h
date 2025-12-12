#pragma once
#include <vector>

namespace hexify {

// Mixed aperture-3/4 hierarchical hexagon grid
// Takes an aperture sequence to handle mixed resolutions

// Quantize to (i,j) at final resolution given aperture sequence
// ap_seq: vector of apertures for each level (3 or 4)
void hex_quantize_ap34(double icosa_triangle_x, double icosa_triangle_y,
                       const std::vector<int>& ap_seq,
                       long long& out_i, long long& out_j);

// Get center of cell (i,j) given aperture sequence
void hex_center_ap34(long long i, long long j,
                     const std::vector<int>& ap_seq,
                     double& out_cx, double& out_cy);

// Get corners of cell (i,j) given aperture sequence
void hex_corners_ap34(long long i, long long j,
                      const std::vector<int>& ap_seq,
                      double hex_radius,
                      double* out_x, double* out_y);  // 6 vertices

} // namespace hexify
