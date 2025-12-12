#pragma once
#include <utility>
#include <string>

namespace hexify {

// Precision presets: "fast", "default", "high", "ultra"
void snyder_inv_set_precision(const std::string& mode = "",
                              double tol_override = -1.0,
                              int    max_iters_override = -1);

std::pair<double,double> snyder_inv_get_precision(); // {tol, max_iters}

void snyder_inv_set_verbose(bool v = true);

// Returns {calls, iters_total, iters_max, capped} and resets stats
std::tuple<int,int,int,int> snyder_inv_get_stats_and_reset();

/**
 * Exact inverse Snyder triangle projection on a fixed face.
 * Input: face-plane (x,y) in normalized coordinates (same as your forward),
 *        triangle/face index 0..19.
 * Optional per-call overrides: non-negative tol/max_iters override globals.
 * Returns (lon_deg, lat_deg).
 */
std::pair<double,double> face_xy_to_ll(double x, double y, int face,
                                       double tol_override = -1.0,
                                       int    max_iters_override = -1);

} // namespace hexify
