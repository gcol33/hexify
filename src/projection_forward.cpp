#include "projection_forward.h"
#include "projection_fuller.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <array>

namespace hexify {

static FaceProjection g_projection = FaceProjection::ISEA;

void use_projection(FaceProjection p) { g_projection = p; }

FaceProjection active_projection() { return g_projection; }

// A point projected to a face, and whether it lies on that face: within DH
// of the face centre and inside the sector boundary.
struct ProjectionWithStatus {
  double x, y;
  bool valid;
};

// The face projection `proj` of a point onto one face. With `validate`, a
// point off the face comes back invalid; without it, the face's formulas are
// applied as they stand.
static ProjectionWithStatus project_core(const Geo& geo, const PolyData& ico_data,
                                         int face, bool validate, FaceProjection proj) {
  const SnyderParams& sp = ico_data.topo->snyder;
  const double glon = geo.lon;
  const double glat = geo.lat;

  const double center_sinlat = ico_data.center_sinlat[face];
  const double center_coslat = ico_data.center_coslat[face];
  const double center_lon    = ico_data.center_lon[face];
  const double face_azimuth = ico_data.face_azimuth_offset[face];

  const double cosLat = std::cos(glat);
  const double sinLat = std::sin(glat);

  double tmp = center_sinlat * sinLat + center_coslat * cosLat * std::cos(glon - center_lon);
  tmp = clampd(tmp, -1.0, 1.0);
  const double z = std::acos(tmp);

  if (validate && z > sp.dh_tolerance) {
    return {0.0, 0.0, false};
  }

  // Point at the face centre: rho = 0 and the azimuth is undefined.
  constexpr double Z_TOLERANCE = 1e-14;
  if (z < Z_TOLERANCE) {
    return {sp.origin_x_off / sp.edge, sp.origin_y_off / sp.edge, true};
  }

  double azimuth = std::atan2(cosLat * std::sin(glon - center_lon),
                              center_coslat * sinLat - center_sinlat * cosLat * std::cos(glon - center_lon))
                   - face_azimuth;
  if (azimuth < 0.0) azimuth += kTwoPi;
  if (azimuth >= kTwoPi) azimuth -= kTwoPi;
  const double face_az = azimuth;

  // Reduce the azimuth to its 120° sector
  const int sector = azimuth_sector(azimuth);
  azimuth -= sector * k2PiOver3;

  const double cos_azimuth = std::cos(azimuth);
  const double sin_azimuth = std::sin(azimuth);

  // Snyder's auxiliary angle for the sector (Snyder notation: δ_z)
  const double dz_angle = std::atan2(sp.tan_el, cos_azimuth + sp.cot_30 * sin_azimuth);

  // The point must lie inside the sector boundary
  if (validate && z > dz_angle + 1e-7) {
    return {0.0, 0.0, false};
  }

  if (proj == FaceProjection::Fuller) {
    const auto xy = fuller_face_xy(z, face_az);
    return {xy.first, xy.second, true};
  }

  // Snyder's angle 'h' - auxiliary spherical angle (Snyder notation: h)
  const double h_arg = clampd(sin_azimuth * sp.sin_g * sp.cos_el - cos_azimuth * sp.cos_g, -1.0, 1.0);
  const double h_angle = std::acos(h_arg);

  // Snyder's accumulated angle from azimuth (Snyder notation: A_G)
  const double AG_angle = azimuth + sp.g_angle + h_angle - kPi;

  // Transformed azimuth in the face plane (Snyder notation: Az')
  double azimuth_transformed = std::atan2(2.0 * AG_angle, sp.r1_squared * sp.tan_el * sp.tan_el - 2.0 * AG_angle * sp.cot_30);

  // Snyder's 'f' scale factor (Snyder notation: f)
  const double denom = 2.0 * (std::cos(azimuth_transformed) + sp.cot_30 * std::sin(azimuth_transformed)) * std::sin(dz_angle / 2.0);
  const double f_scale = (std::fabs(denom) < 1e-15) ? 0.0 : sp.tan_el / denom;

  // Radial distance in face plane (Snyder notation: ρ)
  const double rho = 2.0 * sp.r1 * f_scale * std::sin(z / 2.0);

  // Restore the sector
  azimuth_transformed += sector * k2PiOver3;

  const double x = (rho * std::sin(azimuth_transformed) + sp.origin_x_off) / sp.edge;
  const double y = (rho * std::cos(azimuth_transformed) + sp.origin_y_off) / sp.edge;

  if (validate && (!std::isfinite(x) || !std::isfinite(y))) {
    return {0.0, 0.0, false};
  }
  return {x, y, true};
}

static ProjectionWithStatus project_to_face_with_validation(const Geo& geo, const PolyData& ico_data, int face) {
  return project_core(geo, ico_data, face, /*validate=*/true, active_projection());
}

std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face,
                                         FaceProjection proj) {
  if (face < 0 || face >= ico_data.n_faces()) {
    throw std::invalid_argument("project_to_face: face out of range for the solid");
  }
  const ProjectionWithStatus r = project_core(geo, ico_data, face, /*validate=*/false, proj);
  return {r.x, r.y};
}

std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face) {
  return project_to_face(geo, ico_data, face, active_projection());
}

// Helper: compute great-circle distance from point to face center
static double face_distance(const Geo& point, const PolyData& ico_data, int face) {
  double tmp = ico_data.center_sinlat[face] * std::sin(point.lat) +
               ico_data.center_coslat[face] * std::cos(point.lat) *
               std::cos(ico_data.center_lon[face] - point.lon);
  tmp = clampd(tmp, -1.0, 1.0);
  return std::acos(tmp);
}

ProjectionResult snyder_forward(double lon_deg, double lat_deg) {
  if (!std::isfinite(lon_deg) || !std::isfinite(lat_deg)) {
    throw std::invalid_argument("snyder_forward: lon_deg/lat_deg must be finite (not NA/NaN/Inf)");
  }
  const PolyData& ico_data = poly();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  const int n_faces = ico_data.n_faces();

  // Compute distances to all face centers and sort
  std::array<std::pair<double, int>, kMaxFaces> face_dists;
  for (int i = 0; i < n_faces; ++i) {
    face_dists[i] = {face_distance(g, ico_data, i), i};
  }
  std::sort(face_dists.begin(), face_dists.begin() + n_faces);

  // Try faces in order of distance until one produces valid projection
  const int n_attempts = std::min(5, n_faces);
  for (int attempt = 0; attempt < n_attempts; ++attempt) {
    int face = face_dists[attempt].second;
    auto result = project_to_face_with_validation(g, ico_data, face);
    if (result.valid) {
      return { face, result.x, result.y };
    }
  }

  // Fallback: use closest face with original projection (shouldn't happen)
  int face = face_dists[0].second;
  auto xy = project_to_face(g, ico_data, face);
  return { face, xy.first, xy.second };
}

std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg) {
  const PolyData& ico_data = poly();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  return project_to_face(g, ico_data, face);
}

// Used by the inverse
double snyder_get_face_azimuth_offset(int face) {
  const PolyData& S = poly();
  if (face < 0 || face >= S.n_faces()) throw std::runtime_error("face out of range for the solid");
  return S.face_azimuth_offset[face];
}


} // namespace hexify
