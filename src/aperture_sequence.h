#pragma once
#include <vector>

namespace hexify {

// Mixed-aperture hierarchical hexagon grid.
//
// ap_seq gives the aperture of every resolution level, so its length is
// resolution + 1. ap_seq[0] names the base grid; ap_seq[1..] are the
// refinement steps. Apertures 3, 4 and 7 may appear in any order.

// Quantize to (i,j) at the final resolution of the sequence
void hex_quantize_mixed(double icosa_triangle_x, double icosa_triangle_y,
                        const std::vector<int>& ap_seq,
                        long long& out_i, long long& out_j);

// Get center of cell (i,j) at the final resolution of the sequence
void hex_center_mixed(long long i, long long j,
                      const std::vector<int>& ap_seq,
                      double& out_cx, double& out_cy);

// Get corners of cell (i,j) at the final resolution of the sequence
void hex_corners_mixed(long long i, long long j,
                       const std::vector<int>& ap_seq,
                       double hex_radius,
                       double* out_x, double* out_y);  // 6 vertices

} // namespace hexify
