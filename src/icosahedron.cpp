#include "icosahedron.h"
#include "constants.h"
#include <array>
#include <cmath>

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
  inline Geo coordtrans(const Geo& newNPold, const Geo& ptold, double lon0) {
    // Precompute trigonometric values for clarity
    const double sin_np_lat = std::sin(newNPold.lat);
    const double cos_np_lat = std::cos(newNPold.lat);
    const double sin_pt_lat = std::sin(ptold.lat);
    const double cos_pt_lat = std::cos(ptold.lat);
    const double dlon = newNPold.lon - ptold.lon;
    const double cos_dlon = std::cos(dlon);

    // Compute new latitude using spherical law of cosines
    // cos(new_lat) = sin(np_lat)*sin(pt_lat) + cos(np_lat)*cos(pt_lat)*cos(dlon)
    double cos_new_lat = sin_np_lat * sin_pt_lat + cos_np_lat * cos_pt_lat * cos_dlon;
    cos_new_lat = clampd(cos_new_lat, -1.0, 1.0);
    const double ptnew_lat = std::acos(cos_new_lat);

    // Handle poles (longitude undefined)
    constexpr double POLE_TOLERANCE = kPrecision * 100000;
    const bool at_pole = (std::abs(ptnew_lat) < POLE_TOLERANCE) ||
                         (std::abs(ptnew_lat - kPi) < POLE_TOLERANCE);

    double ptnew_lon = 0.0;
    if (!at_pole) {
      // Compute new longitude using spherical trig identity
      const double sin_new_lat = std::sin(ptnew_lat);
      double cos_new_lon = (sin_pt_lat * cos_np_lat - cos_pt_lat * sin_np_lat * cos_dlon) / sin_new_lat;
      cos_new_lon = clampd(cos_new_lon, -1.0, 1.0);

      ptnew_lon = std::acos(cos_new_lon);

      // Determine longitude sign based on original position
      const double lon_diff = ptold.lon - newNPold.lon;
      if (0.0 <= lon_diff && lon_diff < kPi) {
        ptnew_lon = -ptnew_lon + lon0;
      } else {
        ptnew_lon = ptnew_lon + lon0;
      }
      ptnew_lon = wrap_lon(ptnew_lon);
    }

    // Return co-latitude converted to latitude
    return Geo(ptnew_lon, kPiOver2 - ptnew_lat);
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
    vertsnew[i]   = Geo(wrap_lon(-S_az + deg2rad(72.0 * (i-1))),     deg2rad(26.565051177));
    vertsnew[i+5] = Geo(wrap_lon(-S_az + deg2rad(36.0 + 72.0*(i-1))), -deg2rad(26.565051177));
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

int which_face(double lon_deg, double lat_deg) {
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
