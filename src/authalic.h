#pragma once
// Geodetic and authalic latitude on an ellipsoid of revolution.
//
// The authalic latitude beta of a geodetic latitude phi is the latitude on
// the sphere of the ellipsoid's area (the authalic sphere) below which that
// sphere holds the same area as the ellipsoid below phi:
//   sin(beta) = q(phi) / q_p,
//   q(phi) = (1 - e^2) [sin(phi) / (1 - e^2 sin^2(phi)) + atanh(e sin(phi)) / e],
// q_p = q(pi / 2), e the eccentricity. With longitude kept, (lon, phi) ->
// (lon, beta) carries the ellipsoid onto its authalic sphere with every area
// kept, so an equal-area grid on the sphere is equal-area on the ellipsoid.
//
// A grid's ellipsoid is made active with its solid (activate_icosa()); with
// none active the sphere is the earth model and both conversions return
// their argument unchanged. Lon/lat enter and leave the projection in
// geodetic latitude through to_sphere_lat() and to_geodetic_lat(); the
// projection and everything below it work on the sphere.

namespace hexify {

// Makes the ellipsoid of flattening f active; f = 0 is the sphere.
void use_ellipsoid(double flattening);

// The flattening of the active ellipsoid, 0 for the sphere.
double active_flattening();

// Geodetic -> authalic and authalic -> geodetic latitude in radians on the
// active ellipsoid; the identity on the sphere.
double to_sphere_lat(double lat_rad);
double to_geodetic_lat(double lat_rad);

// The same in degrees.
double to_sphere_lat_deg(double lat_deg);
double to_geodetic_lat_deg(double lat_deg);

// The conversions in closed form on the ellipsoid of flattening f: the
// forward from q(phi), the inverse by Newton's method on it. The active
// ellipsoid's series is fitted to these; they are kept to check it.
double authalic_lat_exact(double phi, double flattening);
double geodetic_lat_exact(double beta, double flattening);

// The authalic radius over the semi-major axis, sqrt(q_p / 2).
double authalic_radius_ratio(double flattening);

// The local scale of the map from the ellipsoid of flattening f onto its
// authalic sphere at geodetic latitude phi: lengths along the parallel are
// multiplied by the returned k and along the meridian by 1 / k.
double authalic_parallel_scale(double phi, double flattening);

} // namespace hexify
