#pragma once
// hex_ap4.h - Aperture-4 hexagon grid

namespace hexify {

// DGGRID aperture 4: Always uses Class I with scale factor 2
// Simpler than aperture 3 (no Class I/II alternation)

// Quantify to (i,j) at a specific resolution
void hex_quantify_ap4(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j);

// Get center of cell (i,j) at resolution
void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);

// Get corners of cell (i,j) at resolution
void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y);  // 6 vertices

} // namespace hexify
