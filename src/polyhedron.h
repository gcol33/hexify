#pragma once
// The base solid a grid is built on, and where it sits on the sphere.
//
// A solid with triangular faces carries a hexagonal grid. Its faces pair into
// diamonds ("quads") across shared edges; each quad has its own (i, j) frame
// and owns the cells whose centres fall in its half-open box. Quad q is the
// solid's vertex q: a diamond quad's origin corner is that vertex, and the two
// vertices no diamond starts at are single-cell quads of their own, quad 0
// and quad n_verts - 1 (the "poles"). The icosahedron has 10 diamonds and the
// octahedron 4. The tetrahedron has no such tiling (an origin vertex is the
// source of exactly two faces, so every vertex of valence 3 must be a pole,
// and the tetrahedron has four), so it carries the face projection only.
//
// Everything that follows from the solid alone -- its faces, Snyder's
// projection constants, the quad tables -- is a SolidTopology, built once per
// solid from its face list. Where the solid sits on the sphere is a PolyData,
// built once per (solid, orientation) and made active by use_solid().

#include <array>
#include <vector>

namespace hexify {

// ---- Basic geographic structure ----
struct Geo {
  double lon; // radians
  double lat; // radians
  Geo() : lon(0.0), lat(0.0) {}
  Geo(double lo, double la) : lon(lo), lat(la) {}
};

// ---- Utility functions ----
double deg2rad(double d);
double rad2deg(double r);
double clampd(double x, double a, double b);
double wrap_lon(double lon_rad);

// ---- Solids ----
enum class Solid { Icosahedron = 0, Octahedron = 1, Tetrahedron = 2 };

constexpr int kMaxFaces = 20;
constexpr int kMaxVerts = 12;

// Where a solid sits on the sphere: vertex 0 and the azimuth of vertex 1 seen
// from it, in degrees (DGGRID's dggs_vert0_lon, dggs_vert0_lat and
// dggs_vert0_azimuth). The defaults are the standard ISEA orientation.
struct Orientation {
  double vert0_lon_deg = 11.25;
  double vert0_lat_deg = 58.282525588538995;
  double azimuth_deg   = 0.0;
};

bool operator==(const Orientation& a, const Orientation& b);

// Snyder's (1992) equal-area projection of one triangular face. g is the
// spherical angle from a face centre to its vertices and G half the angle at
// a vertex, 180 / valence degrees; theta, half the plane angle at the face
// centre, is 30 degrees on every triangular face. The plane triangle has unit
// edge after scaling by 'edge', whose triangle area is the face's share of the
// unit sphere, 4 pi / n_faces.
struct SnyderParams {
  double el_angle;       // g
  double g_angle;        // G
  double tan_el, cos_el, sin_g, cos_g, cot_30;
  double edge;
  double origin_x_off;   // half the edge
  double origin_y_off;   // the triangle's inradius
  double r1, r1_squared; // Snyder's R' and its square
  double dh_tolerance;   // largest distance from a face centre, with slack
};

// The face a region of the quad plane around a quad's origin lies on, and how
// a quad-plane point there maps to that face's triangle coordinates:
// rotate (point + trans) by -60 * rot60 degrees. The six regions are the
// 60-degree wedges at the origin (see compute_subtriangle()); a region the
// solid's angular deficit leaves empty has keep = false.
struct QuadRegion {
  int face;
  double trans_x, trans_y;
  int rot60;
  bool keep;
};

// How a face sits in its quad: rotate its triangle coordinates by
// 60 * rotations degrees counter-clockwise, then subtract the offset.
struct FacePlacement {
  int quad;
  int rotations;
  double offset_x, offset_y;
};

// The map carrying a quad coordinate across one edge of its quad's box into
// the quad on the other side. With d the coordinate's signed distance past the
// edge along the crossed axis (i - top for the far i edge, i for the near one)
// and 'along' the other coordinate, each new coordinate is
// k[0] * top + k[1] * along + k[2] * d. 'pole' names the single-cell quad of
// the corner where a far edge starts (along = 0) when that corner is a pole,
// and is -1 otherwise.
struct QuadEdgeMap {
  int quad;
  int pole;
  int k[2][3];
};

// A face's place in the PLANE layout of the unfolded solid: rotate by
// 60 * rot60 degrees, then add the offset.
struct PlaneTriLayout {
  int rot60;
  double offset_x;
  double offset_y;
};

// A vertex in the solid's standard frame, vertex 0 at the north pole: its
// longitude and latitude in degrees, and whether the orientation's azimuth
// turns it (a pole of the frame has no longitude to turn).
struct StdVertex {
  double lon_deg;
  double lat_deg;
  bool turns;
};

// Quad edges in the order QuadEdgeMap tables keep them: the near i edge
// (i < 0), the near j edge (j < 0), the far i edge and the far j edge.
enum QuadEdge { kEdgeLeft = 0, kEdgeDown = 1, kEdgeRight = 2, kEdgeUp = 3 };

// Corners of a quad's box in the order SolidTopology::corner keeps them:
// (0, 0), (1, 0), (1, 1), (0, 1) in units of the quad's edge.
enum QuadCorner { kCornerOrigin = 0, kCornerI = 1, kCornerFar = 2, kCornerJ = 3 };

struct SolidTopology {
  Solid solid;
  const char* name;
  int n_faces;
  int n_verts;          // one quad per vertex
  bool has_quads;       // false where no diamond tiling exists
  std::array<std::array<int, 3>, kMaxFaces> faces;  // counter-clockwise from outside
  std::array<int, kMaxVerts> valence;
  std::array<std::vector<int>, kMaxVerts> neighbors;  // the vertex graph, sorted
  SnyderParams snyder;
  std::array<StdVertex, kMaxVerts> std_verts;
  Orientation default_orientation;
  std::array<PlaneTriLayout, kMaxFaces> plane;
  // The quad scheme; read only when has_quads.
  std::array<FacePlacement, kMaxFaces> placement;
  std::array<std::array<int, 4>, kMaxVerts> corner;           // QuadCorner order
  std::array<std::array<QuadRegion, 6>, kMaxVerts> region;
  std::array<std::array<QuadEdgeMap, 4>, kMaxVerts> edge;     // QuadEdge order

  int n_quads() const { return n_verts; }
  int n_diamonds() const { return has_quads ? n_verts - 2 : 0; }
  int south_pole() const { return n_verts - 1; }
  bool is_pole(int quad) const { return quad == 0 || quad == n_verts - 1; }
};

// The topology of a solid, built on first use
const SolidTopology& solid_topology(Solid s);

// ---- A solid placed on the sphere ----
struct PolyData {
  const SolidTopology* topo;
  std::array<Geo, kMaxFaces> centers;              // face centres (radians)
  std::array<double, kMaxFaces> center_sinlat;
  std::array<double, kMaxFaces> center_coslat;
  std::array<double, kMaxFaces> center_lon;
  std::array<double, kMaxFaces> face_azimuth_offset;  // towards each face's vertex 0
  std::array<Geo, kMaxVerts> verts;                 // vertices (radians)
  // The flat face as an affine map of its triangle coordinates: the point
  // (tx, ty) of face f sits at solid_origin[f] + tx * solid_x[f] + ty * solid_y[f]
  // on the solid inscribed in the unit sphere.
  std::array<std::array<double, 3>, kMaxFaces> solid_origin;
  std::array<std::array<double, 3>, kMaxFaces> solid_x;
  std::array<std::array<double, 3>, kMaxFaces> solid_y;
  // Each face of the spherical solid is the region on the inner side of the
  // three great circles through its edges: n . p >= 0 for the three normals.
  std::array<std::array<std::array<double, 3>, 3>, kMaxFaces> edge_normal;

  int n_faces() const { return topo->n_faces; }
  int n_verts() const { return topo->n_verts; }
};

// ---- Activation ----

// Sets the default orientation of the icosahedron -- the one read by calls
// that name none -- and makes it the active solid.
void build_icosa_full(double vert0_lon_deg = 11.25,
                      double vert0_lat_deg = 58.282525588538995,
                      double azimuth_deg   = 0.0);

// Makes solid `s` placed by `o` the one poly() and every query below read.
// Each table is built once and kept, so switching between grids costs a
// lookup.
void use_solid(Solid s, const Orientation& o);

// Makes solid `s` in its default orientation the active one.
void use_default_solid(Solid s);

// Makes the default icosahedron the active solid.
void use_default_orientation();

// The orientation calls that name none read for solid `s`
Orientation default_orientation(Solid s);

// ---- Accessors ----

// The active solid on the sphere
const PolyData& poly();

// The active solid's topology
inline const SolidTopology& topo() { return *poly().topo; }

// Face centres of the active solid, in radians
const std::array<Geo, kMaxFaces>& face_centers();

// Per-face azimuth offset (radians)
double get_face_azimuth_offset(int face);

// The face a point (lon_deg, lat_deg) lies on
int which_face(double lon_deg, double lat_deg);

// The point (tx, ty) of a face on the flat solid, as xyz
void face_tri_to_solid(int face, double tx, double ty, double out[3]);

// The point (tx, ty) of a face on the unit sphere, as xyz
void face_tri_to_sphere(int face, double tx, double ty, double out[3]);

// The point (tx, ty) of a face in the PLANE layout of the unfolded solid
void face_tri_to_plane(int face, double tx, double ty, double& px, double& py);

} // namespace hexify
