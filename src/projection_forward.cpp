#include "projection_forward.h"
#include "projection_fuller.h"
#include "projection_ivea.h"
#include "authalic.h"
#include "constants.h"
#include "dual.h"
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

// A point's azimuth within its 120-degree sector of a face: the sector, the
// azimuth from the sector's first vertex with its cosine and sine, and
// Snyder's dz, the arc from the face centre to the face edge along it.
template <class T>
struct SectorAzimuth {
  int sector;
  T azimuth, cos_az, sin_az, dz;
};

// The sector of azimuth az, from the face's first vertex in [0, 2 pi)
template <class T>
static inline SectorAzimuth<T> sector_azimuth(const SnyderParams& sp, T az) {
  using std::atan2;
  using std::cos;
  using std::sin;
  SectorAzimuth<T> s;
  s.sector = azimuth_sector(value_of(az));
  s.azimuth = az - s.sector * k2PiOver3;
  s.cos_az = cos(s.azimuth);
  s.sin_az = sin(s.azimuth);
  s.dz = atan2(sp.tan_el, s.cos_az + sp.cot_30 * s.sin_az);
  return s;
}

// Snyder's (1992) adjustment of Lambert's azimuthal equal-area projection
// about the face centre (eqs. 6-11): the point at azimuth Az within its
// sector goes to the plane azimuth Az' within the sector, which cuts off the
// same area on the plane triangle as Az on the spherical one, and its radius
// 2 sin(z / 2) is multiplied by f, which puts the face edge on the plane
// triangle's edge. f does not depend on z.
template <class T>
static inline void snyder_adjust(const SnyderParams& sp, const SectorAzimuth<T>& s,
                                 T& az_prime, T& f) {
  using std::acos;
  using std::atan2;
  using std::cos;
  using std::sin;
  // Snyder's spherical angle H at the point where the ray at Az meets the edge
  const T h_arg = clamp_to(s.sin_az * sp.sin_g * sp.cos_el - s.cos_az * sp.cos_g, -1.0, 1.0);
  const T h_angle = acos(h_arg);
  // The spherical excess A_G of the triangle cut off by the ray
  const T AG_angle = s.azimuth + sp.g_angle + h_angle - kPi;
  az_prime = atan2(2.0 * AG_angle, sp.r1_squared * sp.tan_el * sp.tan_el - 2.0 * AG_angle * sp.cot_30);
  const T denom = 2.0 * (cos(az_prime) + sp.cot_30 * sin(az_prime)) * sin(s.dz / 2.0);
  f = (std::fabs(value_of(denom)) < 1e-15) ? T(0.0) : sp.tan_el / denom;
}

// Snyder's point on the face plane, (u, v) from the face centre in units of
// the unit sphere along the plane triangle's x and y axes: radius
// rho = 2 R' f sin(z / 2) at the plane azimuth az_prime within `sector`.
template <class T>
static inline void snyder_plane(const SnyderParams& sp, T z, int sector, T az_prime, T f,
                                T& u, T& v) {
  using std::cos;
  using std::sin;
  const T rho = 2.0 * sp.r1 * f * sin(z / 2.0);
  const T a = az_prime + sector * k2PiOver3;
  u = rho * sin(a);
  v = rho * cos(a);
}

// The face projection `proj` of the point at arc z from a face centre and
// azimuth az from the face's first vertex, both in radians, az in [0, 2 pi),
// to face-plane (x, y). With `validate`, a point past the sector boundary
// returns false. T is double, or Dual2 to carry derivatives in z and az.
template <class T>
static bool face_xy_from_polar(const SolidTopology& t, T z, T az, bool validate,
                               FaceProjection proj, T& x, T& y) {
  const SnyderParams& sp = t.snyder;
  const SectorAzimuth<T> s = sector_azimuth(sp, az);

  // The point must lie inside the sector boundary
  if (validate && value_of(z) > value_of(s.dz) + 1e-7) return false;

  if (proj == FaceProjection::Fuller) {
    const auto xy = fuller_face_xy(z, az);
    x = xy.first;
    y = xy.second;
    return true;
  }
  if (proj == FaceProjection::IVEA) {
    const auto xy = ivea_face_xy(t.vgc, z, az);
    x = xy.first;
    y = xy.second;
    return true;
  }

  T az_prime, f, u, v;
  snyder_adjust(sp, s, az_prime, f);
  snyder_plane(sp, z, s.sector, az_prime, f, u, v);
  x = (u + sp.origin_x_off) / sp.edge;
  y = (v + sp.origin_y_off) / sp.edge;
  return true;
}

// The point at arc z and azimuth az from a face centre at step `stage` of
// Lambert's construction, as (u, v) on the plane of that step in units of the
// unit sphere, along the plane triangle's x and y axes from the face centre.
// ConstructionStage::Face is the face projection `proj`.
template <class T>
static void stage_uv(const SolidTopology& t, T z, T az, ConstructionStage stage,
                     FaceProjection proj, T& u, T& v) {
  using std::cos;
  using std::sin;
  const SnyderParams& sp = t.snyder;
  if (stage == ConstructionStage::Lambert) {
    const T rho = 2.0 * sin(z / 2.0);
    u = rho * sin(az);
    v = rho * cos(az);
    return;
  }
  if (stage == ConstructionStage::Nudge) {
    const SectorAzimuth<T> s = sector_azimuth(sp, az);
    T az_prime, f;
    snyder_adjust(sp, s, az_prime, f);
    snyder_plane(sp, z, s.sector, az_prime, f, u, v);
    u = u / sp.r1;
    v = v / sp.r1;
    return;
  }
  T x, y;
  face_xy_from_polar(t, z, az, /*validate=*/false, proj, x, y);
  u = x * sp.edge - sp.origin_x_off;
  v = y * sp.edge - sp.origin_y_off;
}

// Arc z from a face centre and azimuth from the face's first vertex, in
// [0, 2 pi), of a point.
static void face_polar(const Geo& geo, const PolyData& ico_data, int face,
                       double& z, double& azimuth) {
  const double glon = geo.lon;
  const double glat = geo.lat;

  const double center_sinlat = ico_data.center_sinlat[face];
  const double center_coslat = ico_data.center_coslat[face];
  const double center_lon    = ico_data.center_lon[face];

  const double cosLat = std::cos(glat);
  const double sinLat = std::sin(glat);

  // The haversine form keeps z precise near the centre, where acos of the
  // dot product loses half the digits.
  const double h_lat = std::sin(0.5 * (glat - ico_data.centers[face].lat));
  const double h_lon = std::sin(0.5 * (glon - center_lon));
  const double hav = clampd(h_lat * h_lat + center_coslat * cosLat * h_lon * h_lon, 0.0, 1.0);
  z = 2.0 * std::asin(std::sqrt(hav));

  azimuth = std::atan2(cosLat * std::sin(glon - center_lon),
                       center_coslat * sinLat - center_sinlat * cosLat * std::cos(glon - center_lon))
            - ico_data.face_azimuth_offset[face];
  if (azimuth < 0.0) azimuth += kTwoPi;
  if (azimuth >= kTwoPi) azimuth -= kTwoPi;
}

// The face projection `proj` of a point onto one face. With `validate`, a
// point off the face comes back invalid; without it, the face's formulas are
// applied as they stand.
static ProjectionWithStatus project_core(const Geo& geo, const PolyData& ico_data,
                                         int face, bool validate, FaceProjection proj) {
  const SnyderParams& sp = ico_data.topo->snyder;
  double z, azimuth;
  face_polar(geo, ico_data, face, z, azimuth);

  if (validate && z > sp.dh_tolerance) {
    return {0.0, 0.0, false};
  }

  // Point at the face centre: rho = 0 and the azimuth is undefined.
  constexpr double Z_TOLERANCE = 1e-14;
  if (z < Z_TOLERANCE) {
    return {sp.origin_x_off / sp.edge, sp.origin_y_off / sp.edge, true};
  }

  double x, y;
  if (!face_xy_from_polar(*ico_data.topo, z, azimuth, validate, proj, x, y)) {
    return {0.0, 0.0, false};
  }
  if (validate && (!std::isfinite(x) || !std::isfinite(y))) {
    return {0.0, 0.0, false};
  }
  return {x, y, true};
}

static void require_construction(ConstructionStage stage) {
  if (stage == ConstructionStage::Nudge && active_projection() != FaceProjection::ISEA) {
    throw std::invalid_argument("Lambert's construction with Snyder's adjustment is "
                                "the ISEA projection's");
  }
}

FaceScale face_scale(const Geo& geo, int face, ConstructionStage stage) {
  const PolyData& ico_data = poly();
  if (face < 0 || face >= ico_data.n_faces()) {
    throw std::invalid_argument("face_scale: face out of range for the solid");
  }
  require_construction(stage);
  double z, azimuth;
  face_polar(geo, ico_data, face, z, azimuth);
  // At the face centre the azimuth is undefined and the projection has a
  // limit along each direction; read it just off the centre.
  constexpr double kNearCentre = 1e-9;
  if (z < kNearCentre) z = kNearCentre;

  Dual2 u, v;
  stage_uv(*ico_data.topo, Dual2(z, 1.0, 0.0), Dual2(azimuth, 0.0, 1.0), stage,
           active_projection(), u, v);
  // Plane lengths in units of the unit sphere: on the face plane the plane
  // triangle then has the face's area.
  const double across = 1.0 / std::sin(z);
  FaceScale out;
  out.j[0][0] = u.d[0];
  out.j[1][0] = v.d[0];
  out.j[0][1] = u.d[1] * across;
  out.j[1][1] = v.d[1] * across;
  return out;
}

ConstructionPoint construction_point(const Geo& geo, int face) {
  const PolyData& ico_data = poly();
  if (face < 0 || face >= ico_data.n_faces()) {
    throw std::invalid_argument("construction_point: face out of range for the solid");
  }
  if (active_projection() != FaceProjection::ISEA) {
    throw std::invalid_argument("Lambert's construction is the ISEA projection's");
  }
  const SolidTopology& t = *ico_data.topo;
  const SnyderParams& sp = t.snyder;
  ConstructionPoint out;
  face_polar(geo, ico_data, face, out.z, out.az);
  stage_uv(t, out.z, out.az, ConstructionStage::Lambert, FaceProjection::ISEA,
           out.lambert_u, out.lambert_v);
  const SectorAzimuth<double> s = sector_azimuth(sp, out.az);
  double az_prime;
  snyder_adjust(sp, s, az_prime, out.f);
  snyder_plane(sp, out.z, s.sector, az_prime, out.f, out.plane_u, out.plane_v);
  out.az_prime = az_prime + s.sector * k2PiOver3;
  if (out.az_prime < 0.0) out.az_prime += kTwoPi;
  if (out.az_prime >= kTwoPi) out.az_prime -= kTwoPi;
  out.tx = (out.plane_u + sp.origin_x_off) / sp.edge;
  out.ty = (out.plane_v + sp.origin_y_off) / sp.edge;
  return out;
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

// The forward projection of a point on the sphere, in radians
static ProjectionResult forward_on_sphere(const Geo& g) {
  const PolyData& ico_data = poly();
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

static void require_finite(double lon_deg, double lat_deg) {
  if (!std::isfinite(lon_deg) || !std::isfinite(lat_deg)) {
    throw std::invalid_argument("snyder_forward: lon_deg/lat_deg must be finite (not NA/NaN/Inf)");
  }
}

ProjectionResult snyder_forward(double lon_deg, double lat_deg) {
  require_finite(lon_deg, lat_deg);
  return forward_on_sphere(Geo(deg2rad(lon_deg), to_sphere_lat(deg2rad(lat_deg))));
}

ProjectionResult snyder_forward_sphere(double lon_deg, double lat_deg) {
  require_finite(lon_deg, lat_deg);
  return forward_on_sphere(Geo(deg2rad(lon_deg), deg2rad(lat_deg)));
}

std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg) {
  const PolyData& ico_data = poly();
  const Geo g(deg2rad(lon_deg), to_sphere_lat(deg2rad(lat_deg)));
  return project_to_face(g, ico_data, face);
}

// Used by the inverse
double snyder_get_face_azimuth_offset(int face) {
  const PolyData& S = poly();
  if (face < 0 || face >= S.n_faces()) throw std::runtime_error("face out of range for the solid");
  return S.face_azimuth_offset[face];
}


} // namespace hexify
