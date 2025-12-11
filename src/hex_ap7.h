#pragma once
// hex_ap7.h - Aperture-7 hierarchical hexagon grid

namespace hexify {

// DGGRID aperture 7: Uses Class III with alternating variants
// Even resolutions (Class III-I): substrate is sqrt(7)× finer, Class I surrogate
// Odd resolutions (Class III-II): substrate is sqrt(21)× finer, Class II surrogate
// All substrates use Class I grid, rotated by M_AP7_ROT_DEG (~19.1°)

// Quantify to (i,j) at a specific resolution
void hex_quantify_ap7(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j);

// Get center of cell (i,j) at resolution
void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);

// Get corners of cell (i,j) at resolution
void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y);  // 6 vertices

} // namespace hexify
