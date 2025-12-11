#include "snyder_forward.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

// TU-local pi (internal linkage; won’t collide)
static const double kPI = std::acos( -1.0 );

// DGGRID Snyder constants (plain const; std::tan/cos/sin aren’t constexpr)
static const double R1    = 0.9103832815;
static const double R1S   = R1 * R1;
static const double DH    = deg2rad(37.37736814); // use core_icosa’s deg2rad
static const double GH    = deg2rad(36.0);
static const double cot30 = 1.0 / std::tan( deg2rad(30.0) );
static const double tanDH = std::tan(DH);
static const double cosDH = std::cos(DH);
static const double sinGH = std::sin(GH);
static const double cosGH = std::cos(GH);

// forward’s screen/face offsets and edge scale
static const double originXOff = 0.6022955029;
static const double originYOff = 0.3477354707;
static const double icosaEdge  = 2.0 * originXOff;

std::pair<double,double> sllxy(const Geo& geo, const IcosaData& sph, int nTri) {
  const double glon = geo.lon;
  const double glat = geo.lat;

  const double cent_sin = sph.cent_sinlat[nTri];
  const double cent_cos = sph.cent_coslat[nTri];
  const double clon     = sph.cent_lon[nTri];
  const double dz0      = sph.dazh[nTri];

  const double cosLat = std::cos(glat);
  const double sinLat = std::sin(glat);

  double tmp = cent_sin * sinLat + cent_cos * cosLat * std::cos(glon - clon);
  tmp = clampd(tmp, -1.0, 1.0); // from core_icosa.h
  const double z = std::acos(tmp);

  double azh = std::atan2(cosLat * std::sin(glon - clon),
                          cent_cos * sinLat - cent_sin * cosLat * std::cos(glon - clon))
               - dz0;

  if (azh < 0.0) azh += 2.0 * kPI;
  const double azh0 = azh;

  if ( (120.0 * kPI / 180.0) <= azh && azh <= (240.0 * kPI / 180.0) ) {
    azh -= (120.0 * kPI / 180.0);
  }
  if (azh > (240.0 * kPI / 180.0)) {
    azh -= (240.0 * kPI / 180.0);
  }

  const double cosAzh = std::cos(azh);
  const double sinAzh = std::sin(azh);
  const double dz = std::atan2(tanDH, cosAzh + cot30 * sinAzh);

  const double h  = std::acos(sinAzh * sinGH * cosDH - cosAzh * cosGH);
  const double ag = azh + GH + h - kPI;
  double azh1 = std::atan2(2.0 * ag, R1S * tanDH * tanDH - 2.0 * ag * cot30);
  const double fh = tanDH / (2.0 * (std::cos(azh1) + cot30 * std::sin(azh1)) * std::sin(dz / 2.0));
  const double ph = 2.0 * R1 * fh * std::sin(z / 2.0);

  if ( (120.0 * kPI / 180.0) <= azh0 && azh0 < (240.0 * kPI / 180.0) ) {
    azh1 += (120.0 * kPI / 180.0);
  }
  if (azh0 >= (240.0 * kPI / 180.0)) {
    azh1 += (240.0 * kPI / 180.0);
  }

  const double x = (ph * std::sin(azh1) + originXOff) / icosaEdge;
  const double y = (ph * std::cos(azh1) + originYOff) / icosaEdge;
  return {x, y};
}

FwdOut snyder_fwd(double lon_deg, double lat_deg) {
  const IcosaData& S = ico();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  const int face = which_face(lon_deg, lat_deg);
  auto xy = sllxy(g, S, face);
  return { face, xy.first, xy.second };
}

std::pair<double,double> snyder_fwd_face(int face, double lon_deg, double lat_deg) {
  const IcosaData& S = ico();
  const Geo g(deg2rad(lon_deg), deg2rad(lat_deg));
  return sllxy(g, S, face);
}

// Used by the inverse
double snyder_face_dazh_rad(int face) {
  if (face < 0 || face >= 20) throw std::runtime_error("face must be 0..19");
  const IcosaData& S = ico();
  return S.dazh[face];
}

} // namespace hexify
