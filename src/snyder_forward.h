#pragma once
#include "core_icosa.h"
#include <utility>

namespace hexify {

// Low-level: face-fixed fwd
std::pair<double,double> sllxy(const Geo& geo, const IcosaData& sph, int nTri);

// High-level
// Returns (face, tx, ty)
struct FwdOut { int face; double tx; double ty; };
FwdOut snyder_fwd(double lon_deg, double lat_deg);

std::pair<double,double> snyder_fwd_face(int face, double lon_deg, double lat_deg);

// NEW: per-face azimuth offset (radians), face = 0..19
double snyder_face_dazh_rad(int face);

} // namespace hexify
