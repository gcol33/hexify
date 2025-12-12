#include "projection_forward.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

// Derived trigonometric values (not constexpr because std::tan/cos/sin aren't constexpr)
static const double COT_30 = 1.0 / std::tan(kPiOver6);  // cot(30°) = 1/tan(π/6)
static const double TAN_EL = std::tan(kSnyderElAngle);
static const double COS_EL = std::cos(kSnyderElAngle);
static const double SIN_G  = std::sin(kSnyderGAngle);
static const double COS_G  = std::cos(kSnyderGAngle);

std::pair<double,double> project_to_face(const Geo& geo, const IcosaData& ico_data, int face) {
  const double glon = geo.lon;
  const double glat = geo.lat;

  const double center_sinlat = ico_data.center_sinlat[face];
  const double center_coslat = ico_data.center_coslat[face];
  const double center_lon    = ico_data.center_lon[face];
  const double face_azimuth = ico_data.face_azimuth_offset[face];

  const double cosLat = std::cos(glat);
  const double sinLat = std::sin(glat);

  double tmp = center_sinlat * sinLat + center_coslat * cosLat * std::cos(glon - center_lon);
  tmp = clampd(tmp, -1.0, 1.0); // from icosahedron.h
  const double z = std::acos(tmp);

  double azimuth = std::atan2(cosLat * std::sin(glon - center_lon),
                              center_coslat * sinLat - center_sinlat * cosLat * std::cos(glon - center_lon))
                   - face_azimuth;

  if (azimuth < 0.0) azimuth += kTwoPi;
  const double azimuth_original = azimuth;

  // Reduce azimuth to [0, 120°) sector
  if (k2PiOver3 <= azimuth && azimuth <= k4PiOver3) {
    azimuth -= k2PiOver3;
  }
  if (azimuth > k4PiOver3) {
    azimuth -= k4PiOver3;
  }

  const double cos_azimuth = std::cos(azimuth);
  const double sin_azimuth = std::sin(azimuth);

  // Snyder's auxiliary angle for the sector (Snyder notation: δ_z)
  const double dz_angle = std::atan2(TAN_EL, cos_azimuth + COT_30 * sin_azimuth);

  // Snyder's angle 'h' - auxiliary spherical angle (Snyder notation: h)
  const double h_angle = std::acos(sin_azimuth * SIN_G * COS_EL - cos_azimuth * COS_G);

  // Snyder's accumulated angle from azimuth (Snyder notation: A_G)
  const double AG_angle = azimuth + kSnyderGAngle + h_angle - kPi;

  // Transformed azimuth in the face plane (Snyder notation: Az')
  double azimuth_transformed = std::atan2(2.0 * AG_angle, kSnyderR1Squared * TAN_EL * TAN_EL - 2.0 * AG_angle * COT_30);

  // Snyder's 'f' scale factor (Snyder notation: f)
  const double f_scale = TAN_EL / (2.0 * (std::cos(azimuth_transformed) + COT_30 * std::sin(azimuth_transformed)) * std::sin(dz_angle / 2.0));

  // Radial distance in face plane (Snyder notation: ρ)
  const double rho = 2.0 * kSnyderR1 * f_scale * std::sin(z / 2.0);

  // Restore to original sector
  if (k2PiOver3 <= azimuth_original && azimuth_original < k4PiOver3) {
    azimuth_transformed += k2PiOver3;
  }
  if (azimuth_original >= k4PiOver3) {
    azimuth_transformed += k4PiOver3;
  }

  const double x = (rho * std::sin(azimuth_transformed) + kSnyderOriginXOff) / kSnyderIcosaEdge;
  const double y = (rho * std::cos(azimuth_transformed) + kSnyderOriginYOff) / kSnyderIcosaEdge;
  return {x, y};
}

ProjectionResult snyder_forward(double lon_deg, double lat_deg) {
  const IcosaData& ico_data = ico();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  const int face = which_face(lon_deg, lat_deg);
  auto xy = project_to_face(g, ico_data, face);
  return { face, xy.first, xy.second };
}

std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg) {
  const IcosaData& ico_data = ico();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  return project_to_face(g, ico_data, face);
}

// Used by the inverse
double snyder_get_face_azimuth_offset(int face) {
  if (face < 0 || face >= 20) throw std::runtime_error("face must be 0..19");
  const IcosaData& S = ico();
  return S.face_azimuth_offset[face];
}

} // namespace hexify
