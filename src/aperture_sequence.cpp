// aperture_sequence.cpp - Mixed-aperture hierarchical hexagon grid
//
// Handles arbitrary sequences of aperture 3, 4 and 7 resolutions, covering
// mixed grids such as ISEA43H (aperture 4 for the coarse levels, aperture 3
// for the fine ones) as well as sequences mixing in aperture 7.
//
// Scale and orientation both come from the sequence: each step multiplies the
// scale by sqrt(aperture) and the lattice generator by that aperture's
// Eisenstein integer. See grid_math.h for the model.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "aperture_sequence.h"
#include "grid_math.h"

namespace hexify {

void hex_quantize_mixed(double icosa_triangle_x, double icosa_triangle_y,
                        const std::vector<int>& ap_seq,
                        long long& out_i, long long& out_j) {
    HexGridForm form = hex_form_sequence(ap_seq);
    quantize_form(form, icosa_triangle_x, icosa_triangle_y, out_i, out_j);
}

void hex_center_mixed(long long i, long long j,
                      const std::vector<int>& ap_seq,
                      double& out_cx, double& out_cy) {
    HexGridForm form = hex_form_sequence(ap_seq);
    center_form(form, i, j, out_cx, out_cy);
}

void hex_corners_mixed(long long i, long long j,
                       const std::vector<int>& ap_seq,
                       double hex_radius,
                       double* out_x, double* out_y) {
    HexGridForm form = hex_form_sequence(ap_seq);
    corners_form(form, i, j, hex_radius, out_x, out_y);
}

} // namespace hexify
