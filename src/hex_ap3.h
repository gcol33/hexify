#pragma once
// hex_ap3.h - Aperture-3 hierarchical hexagon grid (Class I & II)

namespace hexify {

// DGGRID-EXACT: Quantify to (i,j) at a specific resolution
// Automatically handles Class I (even res) and Class II (odd res)
// No digit sequences - just absolute grid coordinates
void hex_quantify_ap3(double tx, double ty, int resolution, 
                      long long& out_i, long long& out_j);

// DGGRID-EXACT: Get center of cell (i,j) at resolution
// Automatically handles Class I (even res) and Class II (odd res)
void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);

// DGGRID-EXACT: Get corners of cell (i,j) at resolution
// Automatically handles Class I (even res) and Class II (odd res)
// Class II hexagons are rotated 30° relative to Class I
void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y);  // 6 vertices

} // namespace hexify
