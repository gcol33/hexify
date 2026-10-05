#pragma once
#include <array>

namespace hexify {

// ---- Basic geographic structure ----
struct Geo {
  double lon; // radians
  double lat; // radians
  Geo() : lon(0.0), lat(0.0) {}
  Geo(double lo, double la) : lon(lo), lat(la) {}
};

// ---- Icosahedron data ----
struct IcosaData {
  std::array<Geo, 20> centers;            // face centers (radians)
  std::array<double, 20> center_sinlat;   // sin(center latitude) per face
  std::array<double, 20> center_coslat;   // cos(center latitude) per face
  std::array<double, 20> center_lon;      // center longitude per face (radians)
  std::array<double, 20> face_azimuth_offset;  // per-face azimuth offsets (radians)
  std::array<Geo, 12> verts;              // icosahedron vertices (radians)
  std::array<std::array<int, 3>, 20> face_verts;  // vertex indices of each face
  // The flat face as an affine map of its triangle coordinates: the point
  // (tx, ty) of face f sits at solid_origin[f] + tx * solid_x[f] + ty * solid_y[f]
  // on the icosahedron inscribed in the unit sphere.
  std::array<std::array<double, 3>, 20> solid_origin;
  std::array<std::array<double, 3>, 20> solid_x;
  std::array<std::array<double, 3>, 20> solid_y;
  bool built = false;
};

// ---- Utility functions ----
double deg2rad(double d);
double rad2deg(double r);
double clampd(double x, double a, double b);
double wrap_lon(double lon_rad);

// ---- Icosahedron construction and queries ----
void build_icosa_full(double vert0_lon_deg = 11.25,
                      double vert0_lat_deg = 58.282525588538995,
                      double azimuth_deg   = 0.0);

// Identify which icosahedral face a point (lon_deg, lat_deg) belongs to
int which_face(double lon_deg, double lat_deg);

// ---- Accessors ----

// Returns the global IcosaData (builds it if needed)
const IcosaData& ico();

// Returns all face centers in radians
const std::array<Geo,20>& face_centers();

// Returns the per-face azimuth offset (in radians)
double get_face_azimuth_offset(int face);

// The point (tx, ty) of a face on the flat icosahedron, as xyz
void face_tri_to_solid(int face, double tx, double ty, double out[3]);

// The point (tx, ty) of a face on the unit sphere, as xyz
void face_tri_to_sphere(int face, double tx, double ty, double out[3]);

// The point (tx, ty) of a face in the PLANE layout of the unfolded icosahedron
void face_tri_to_plane(int face, double tx, double ty, double& px, double& py);

} // namespace hexify
