#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "constants.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {
  constexpr double kPrecision = 1e-15;

  IcosaData g_ico;

  struct Vec3 { double x, y, z; };

  inline Vec3 ll2xyz(const Geo& g) {
    const double cl = std::cos(g.lat);
    return { cl * std::cos(g.lon), cl * std::sin(g.lon), std::sin(g.lat) };
  }

  inline Geo xyz2ll(const Vec3& v_in) {
    const double n = std::sqrt(v_in.x*v_in.x + v_in.y*v_in.y + v_in.z*v_in.z);
    const double x = v_in.x / n, y = v_in.y / n, z = v_in.z / n;
    return Geo(std::atan2(y, x), std::asin(z));
  }

  inline Geo sph_tricen(const Geo tri[3]) {
    const Vec3 a = ll2xyz(tri[0]);
    const Vec3 b = ll2xyz(tri[1]);
    const Vec3 c = ll2xyz(tri[2]);
    const Vec3 v{ a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z };
    return xyz2ll(v);
  }

  /**
   * Great-circle distance between two points on sphere.
   * Uses spherical law of cosines: cos(c) = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
   */
  inline double gc_dist(const Geo& A, const Geo& B) {
    const double sin_lat_A = std::sin(A.lat);
    const double sin_lat_B = std::sin(B.lat);
    const double cos_lat_A = std::cos(A.lat);
    const double cos_lat_B = std::cos(B.lat);
    const double dlon = A.lon - B.lon;

    // Spherical law of cosines
    double cos_dist = sin_lat_A * sin_lat_B + cos_lat_A * cos_lat_B * std::cos(dlon);
    cos_dist = clampd(cos_dist, -1.0, 1.0);

    double dist = std::acos(cos_dist);
    if (dist > kPi) dist = kTwoPi - dist;

    return dist;
  }

  /**
   * Coordinate transformation: rotates point ptold to new coordinate system.
   *
   * newNPold: location of new North Pole in old coordinate system
   * ptold: point to transform (in old coordinates)
   * lon0: longitude offset for new system
   *
   * Mathematical basis: spherical coordinate rotation via great-circle geometry
   */
  //
  // The new co-latitude is the angle between the point and the new pole, and
  // the new longitude is lon0 less the point's azimuth seen from the new pole,
  // measured from the direction of the old north pole and positive eastward.
  // Both come from atan2 on unit vectors, which keeps full precision at every
  // angle.
  inline Geo coordtrans(const Geo& newNPold, const Geo& ptold, double lon0) {
    const Vec3 p = ll2xyz(ptold);
    const Vec3 n = ll2xyz(newNPold);
    // Unit vectors at the new pole: towards the old north pole, and east.
    const double sl = std::sin(newNPold.lat), cl = std::cos(newNPold.lat);
    const double so = std::sin(newNPold.lon), co = std::cos(newNPold.lon);
    const Vec3 north{ -sl * co, -sl * so, cl };
    const Vec3 east{ -so, co, 0.0 };

    const double along = p.x * n.x + p.y * n.y + p.z * n.z;
    const double cx = p.y * n.z - p.z * n.y;
    const double cy = p.z * n.x - p.x * n.z;
    const double cz = p.x * n.y - p.y * n.x;
    const double across = std::sqrt(cx * cx + cy * cy + cz * cz);
    const double colat = std::atan2(across, along);

    // Longitude is undefined at the new pole and its antipode.
    constexpr double POLE_TOLERANCE = kPrecision * 100000;
    double lon = 0.0;
    if (across >= POLE_TOLERANCE) {
      const double az = std::atan2(p.x * east.x + p.y * east.y + p.z * east.z,
                                   p.x * north.x + p.y * north.y + p.z * north.z);
      lon = wrap_lon(lon0 - az);
    }
    return Geo(lon, kPiOver2 - colat);
  }
  
} // anon

// ---- exported helpers (single definitions; others should call these) ----
double deg2rad(double degrees) { return degrees * kDegToRad; }
double rad2deg(double radians) { return radians * kRadToDeg; }
double clampd(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
double wrap_lon(double lon_rad) {
  if (lon_rad >  kPi) lon_rad -= kTwoPi;
  if (lon_rad < -kPi) lon_rad += kTwoPi;
  return lon_rad;
}

// ---- build + queries ----
void build_icosa_full(double vert0_lon_deg, double vert0_lat_deg, double azimuth_deg) {
  const Geo S_pt(deg2rad(vert0_lon_deg), deg2rad(vert0_lat_deg));
  const double S_az = deg2rad(azimuth_deg);

  std::array<Geo,12> vertsnew;
  const Geo newnpold(0.0, S_pt.lat);

  for (int i = 1; i <= 5; ++i) {
    vertsnew[i]   = Geo(wrap_lon(-S_az + deg2rad(72.0 * (i-1))),     deg2rad(kIcosaVertexLatDeg));
    vertsnew[i+5] = Geo(wrap_lon(-S_az + deg2rad(36.0 + 72.0*(i-1))), -deg2rad(kIcosaVertexLatDeg));
  }
  vertsnew[11] = Geo(0.0, -deg2rad(90.0));

  std::array<Geo,12> icoverts;
  icoverts[0] = S_pt;
  for (int i = 1; i < 12; ++i) {
    icoverts[i] = coordtrans(newnpold, vertsnew[i], S_pt.lon);
  }

  static const int faces[20][3] = {
    {0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},
    {6,2,1},{7,3,2},{8,4,3},{9,5,4},{10,1,5},
    {2,6,7},{3,7,8},{4,8,9},{5,9,10},{1,10,6},
    {11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10}
  };

  for (int i = 0; i < 20; ++i) {
    Geo tri[3] = { icoverts[faces[i][0]], icoverts[faces[i][1]], icoverts[faces[i][2]] };
    Geo c = sph_tricen(tri);
    g_ico.centers[i]       = c;
    g_ico.center_sinlat[i] = std::sin(c.lat);
    g_ico.center_coslat[i] = std::cos(c.lat);
    g_ico.center_lon[i]    = c.lon;
  }

  for (int i = 0; i < 20; ++i) {
    const Geo& c  = g_ico.centers[i];
    const Geo& t0 = icoverts[faces[i][0]];
    const double num = std::cos(t0.lat) * std::sin(t0.lon - c.lon);
    const double den = g_ico.center_coslat[i] * std::sin(t0.lat)
                     - std::sin(c.lat) * std::cos(t0.lat) * std::cos(t0.lon - c.lon);
    g_ico.face_azimuth_offset[i] = std::atan2(num, den);
  }

  for (int v = 0; v < 12; ++v) g_ico.verts[v] = icoverts[v];
  for (int i = 0; i < 20; ++i) {
    for (int k = 0; k < 3; ++k) g_ico.face_verts[i][k] = faces[i][k];
  }

  // A face's corners project to its triangle's corners, so the three pairs
  // (triangle coordinates, vertex position) fix the face's affine map.
  for (int i = 0; i < 20; ++i) {
    Vec3 P[3];
    double t[3][2];
    for (int k = 0; k < 3; ++k) {
      P[k] = ll2xyz(icoverts[faces[i][k]]);
      const auto xy = project_to_face(icoverts[faces[i][k]], g_ico, i);
      t[k][0] = xy.first;
      t[k][1] = xy.second;
    }
    const double a = t[1][0] - t[0][0], b = t[2][0] - t[0][0];
    const double c = t[1][1] - t[0][1], d = t[2][1] - t[0][1];
    const double det = a * d - b * c;
    const double inv[2][2] = { {  d / det, -b / det },
                               { -c / det,  a / det } };
    const double E1[3] = { P[1].x - P[0].x, P[1].y - P[0].y, P[1].z - P[0].z };
    const double E2[3] = { P[2].x - P[0].x, P[2].y - P[0].y, P[2].z - P[0].z };
    const double P0[3] = { P[0].x, P[0].y, P[0].z };
    for (int r = 0; r < 3; ++r) {
      g_ico.solid_x[i][r] = E1[r] * inv[0][0] + E2[r] * inv[1][0];
      g_ico.solid_y[i][r] = E1[r] * inv[0][1] + E2[r] * inv[1][1];
      g_ico.solid_origin[i][r] = P0[r] - g_ico.solid_x[i][r] * t[0][0]
                                       - g_ico.solid_y[i][r] * t[0][1];
    }
  }

  g_ico.built = true;
}

const IcosaData& ico() {
  if (!g_ico.built) build_icosa_full();
  return g_ico;
}

const std::array<Geo,20>& face_centers() { return ico().centers; }

// accessor needed by the inverse
double get_face_azimuth_offset(int face) {
  const IcosaData& icosa = ico();
  if (face < 0 || face >= 20) return 0.0;
  return icosa.face_azimuth_offset[face];
}

void face_tri_to_solid(int face, double tx, double ty, double out[3]) {
  const IcosaData& S = ico();
  for (int r = 0; r < 3; ++r) {
    out[r] = S.solid_origin[face][r] + tx * S.solid_x[face][r] + ty * S.solid_y[face][r];
  }
}

void face_tri_to_sphere(int face, double tx, double ty, double out[3]) {
  const auto ll = face_xy_to_ll(tx, ty, face);
  const Vec3 v = ll2xyz(Geo(deg2rad(ll.first), deg2rad(ll.second)));
  out[0] = v.x;
  out[1] = v.y;
  out[2] = v.z;
}

void face_tri_to_plane(int face, double tx, double ty, double& px, double& py) {
  const PlaneTriLayout& layout = kPlaneLayout[face];
  if (layout.rot60 != 0) {
    const double a = layout.rot60 * 60.0 * kDegToRad;
    const double c = std::cos(a), s = std::sin(a);
    const double x = tx * c - ty * s;
    ty = tx * s + ty * c;
    tx = x;
  }
  px = tx + layout.offset_x;
  py = ty + layout.offset_y;
}

int which_face(double lon_deg, double lat_deg) {
  if (!std::isfinite(lon_deg) || !std::isfinite(lat_deg)) {
    throw std::invalid_argument("which_face: lon_deg/lat_deg must be finite (not NA/NaN/Inf)");
  }
  const IcosaData& icosa = ico();
  const Geo point(deg2rad(lon_deg), deg2rad(lat_deg));
  int best = 0;
  double bestd = std::acos(clampd(std::sin(icosa.centers[0].lat)*std::sin(point.lat)
                       + std::cos(icosa.centers[0].lat)*std::cos(point.lat)*std::cos(icosa.centers[0].lon - point.lon),
                       -1.0, 1.0));
  for (int i = 1; i < 20; ++i) {
    const auto& c = icosa.centers[i];
    const double cc = clampd(std::sin(c.lat)*std::sin(point.lat)
                     + std::cos(c.lat)*std::cos(point.lat)*std::cos(c.lon - point.lon), -1.0, 1.0);
    const double d = std::acos(cc);
    if (d < bestd) { best = i; bestd = d; }
  }
  
  return best;
}

} // namespace hexify
