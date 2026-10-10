#include "polyhedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "authalic.h"
#include "constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace hexify {

namespace {
  constexpr double kPrecision = 1e-15;

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

  // Coordinate transformation: rotates point ptold to the frame whose north
  // pole is newNPold. The new co-latitude is the angle between the point and
  // the new pole, and the new longitude is lon0 less the point's azimuth seen
  // from the new pole, measured from the direction of the old north pole and
  // positive eastward. Both come from atan2 on unit vectors, which keeps full
  // precision at every angle.
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

bool operator==(const Orientation& a, const Orientation& b) {
  return a.vert0_lon_deg == b.vert0_lon_deg &&
         a.vert0_lat_deg == b.vert0_lat_deg &&
         a.azimuth_deg == b.azimuth_deg;
}

// ============================================================================
// Solids
// ============================================================================

namespace {

// A solid as data: its faces counter-clockwise seen from outside, how they
// pair into diamonds, and the few choices the quad frames leave open.
struct SolidSpec {
  Solid solid;
  const char* name;
  int n_verts;
  std::vector<std::array<int, 3>> faces;
  // Per diamond quad 1 .. n_verts - 2: the faces in its upper (primary) and
  // lower (secondary) triangle. Empty for a solid without a diamond tiling.
  std::vector<std::array<int, 2>> diamonds;
  // Per quad: how many faces around its origin are laid counter-clockwise
  // from the region at 0-60 degrees before the solid's angular deficit; the
  // rest are laid clockwise from the region at 300-360 degrees.
  std::vector<int> n_ccw;
  // The face at 0-60 degrees in the frame of each pole quad.
  std::array<int, 2> pole_first_face;
  std::vector<StdVertex> std_verts;
  Orientation default_orientation;
  std::vector<PlaneTriLayout> plane;
  double el_angle;   // Snyder's g, from its closed form
  double g_deg;      // Snyder's G
};

SnyderParams snyder_params(double el_angle, double g_deg, int n_faces) {
  SnyderParams s;
  s.el_angle = el_angle;
  s.g_angle = g_deg * kDegToRad;
  s.tan_el = std::tan(el_angle);
  s.cos_el = std::cos(el_angle);
  s.sin_g = std::sin(s.g_angle);
  s.cos_g = std::cos(s.g_angle);
  s.cot_30 = 1.0 / std::tan(kPiOver6);
  // The plane triangle's area (sqrt(3)/4) e^2 is the face's share of the unit
  // sphere, 4 pi / n_faces.
  s.edge = std::sqrt(4.0 * kPi / ((n_faces / 4.0) * std::sqrt(3.0)));
  s.origin_x_off = 0.5 * s.edge;
  s.origin_y_off = s.edge / (2.0 * std::sqrt(3.0));
  // A face vertex lies at plane distance R' tan g from the centre, which is
  // the triangle's circumradius e / sqrt(3).
  s.r1 = s.edge / (std::sqrt(3.0) * std::tan(el_angle));
  s.r1_squared = s.r1 * s.r1;
  s.dh_tolerance = el_angle + 1e-10;
  const double sg = std::sin(el_angle), cg = std::cos(el_angle);
  const double centre[3] = {0.0, 0.0, 1.0};
  const double vert_a[3] = {0.0, sg, cg};
  const double vert_b[3] = {0.5 * std::sqrt(3.0) * sg, -0.5 * sg, cg};
  s.sector = make_snyder_triangle(centre, vert_a, vert_b);
  return s;
}

VertexGcParams vertex_gc_params(double el_angle, double g_deg) {
  VertexGcParams p;
  p.beta = g_deg * kDegToRad;
  p.gamma = kPiOver3;
  p.excess = p.beta + p.gamma - kPiOver2;
  // Right angle at A: cos(beta) = tan(AB) / tan(BC)
  p.tan_ab = std::cos(p.beta) * std::tan(el_angle);
  p.cos_ab = 1.0 / std::sqrt(1.0 + p.tan_ab * p.tan_ab);
  p.bc = el_angle;
  p.sin_bc = std::sin(el_angle);
  p.cos_bc = std::cos(el_angle);
  return p;
}

// The standard ISEA icosahedron: vertex 0 at the frame's north pole, five
// vertices at latitude atan(1/2), five at -atan(1/2) half a step round, and
// vertex 11 at the south pole.
SolidSpec icosahedron_spec() {
  SolidSpec s;
  s.solid = Solid::Icosahedron;
  s.name = "icosahedron";
  s.n_verts = 12;
  s.faces = {
    {0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},
    {6,2,1},{7,3,2},{8,4,3},{9,5,4},{10,1,5},
    {2,6,7},{3,7,8},{4,8,9},{5,9,10},{1,10,6},
    {11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10}
  };
  s.diamonds = {{0,5},{1,6},{2,7},{3,8},{4,9},{10,15},{11,16},{12,17},{13,18},{14,19}};
  s.n_ccw = {4, 4,4,4,4,4, 3,3,3,3,3, 3};
  s.pole_first_face = {0, 18};
  s.std_verts.assign(12, StdVertex{0.0, 90.0, false});
  for (int i = 1; i <= 5; ++i) {
    s.std_verts[i]     = StdVertex{72.0 * (i - 1), kIcosaVertexLatDeg, true};
    s.std_verts[i + 5] = StdVertex{36.0 + 72.0 * (i - 1), -kIcosaVertexLatDeg, true};
  }
  s.std_verts[11] = StdVertex{0.0, -90.0, false};
  s.default_orientation = Orientation{};
  // Five strips of four triangles: the cap and band triangles of an upper
  // quad, then those of a lower quad, 5.5 by 1.73 units in all.
  s.plane = {
    {0, 0.0, 2.0 * kSin60}, {0, 1.0, 2.0 * kSin60}, {0, 2.0, 2.0 * kSin60},
    {0, 3.0, 2.0 * kSin60}, {0, 4.0, 2.0 * kSin60},
    {3, 1.0, 2.0 * kSin60}, {3, 2.0, 2.0 * kSin60}, {3, 3.0, 2.0 * kSin60},
    {3, 4.0, 2.0 * kSin60}, {3, 5.0, 2.0 * kSin60},
    {0, 0.5, kSin60}, {0, 1.5, kSin60}, {0, 2.5, kSin60}, {0, 3.5, kSin60},
    {0, 4.5, kSin60},
    {3, 1.5, kSin60}, {3, 2.5, kSin60}, {3, 3.5, kSin60}, {3, 4.5, kSin60},
    {3, 5.5, kSin60}
  };
  // tan g = 3 - sqrt(5)
  s.el_angle = std::atan(3.0 - std::sqrt(5.0));
  s.g_deg = 36.0;
  return s;
}

// The octahedron: vertex 0 at the frame's north pole, four vertices on the
// equator, vertex 5 at the south pole. Each diamond is a northern face and
// the southern face across its equatorial edge.
SolidSpec octahedron_spec() {
  SolidSpec s;
  s.solid = Solid::Octahedron;
  s.name = "octahedron";
  s.n_verts = 6;
  s.faces = {
    {0,1,2},{0,2,3},{0,3,4},{0,4,1},
    {5,2,1},{5,3,2},{5,4,3},{5,1,4}
  };
  s.diamonds = {{0,4},{1,5},{2,6},{3,7}};
  s.n_ccw = {3, 3,3,3,3, 2};
  s.pole_first_face = {0, 4};
  s.std_verts.assign(6, StdVertex{0.0, 90.0, false});
  for (int i = 1; i <= 4; ++i) s.std_verts[i] = StdVertex{90.0 * (i - 1), 0.0, true};
  s.std_verts[5] = StdVertex{0.0, -90.0, false};
  // Vertex 1 on the prime meridian.
  s.default_orientation = Orientation{0.0, 90.0, 180.0};
  // Four diamonds side by side, each a northern triangle over its southern
  // one, 4 by 1.73 units.
  s.plane = {
    {0, 0.0, kSin60}, {0, 1.0, kSin60}, {0, 2.0, kSin60}, {0, 3.0, kSin60},
    {3, 1.0, kSin60}, {3, 2.0, kSin60}, {3, 3.0, kSin60}, {3, 4.0, kSin60}
  };
  // tan g = sqrt(2)
  s.el_angle = std::atan(std::sqrt(2.0));
  s.g_deg = 45.0;
  return s;
}

// The tetrahedron: vertex 0 at the frame's north pole and three vertices at
// latitude -asin(1/3). Its two diamonds start at vertices 1 and 2 and have
// the opposite edges 1-0 and 2-3 as diagonals; vertices 0 and 3 are the
// poles, each the far corner of one diamond and the side corner of the other.
SolidSpec tetrahedron_spec() {
  SolidSpec s;
  s.solid = Solid::Tetrahedron;
  s.name = "tetrahedron";
  s.n_verts = 4;
  s.faces = {{0,1,2},{0,2,3},{0,3,1},{1,3,2}};
  s.diamonds = {{2,0},{1,3}};
  s.n_ccw = {3, 3,3, 3};
  s.pole_first_face = {0, 3};
  const double lat = -std::asin(1.0 / 3.0) * kRadToDeg;
  s.std_verts.assign(4, StdVertex{0.0, 90.0, false});
  for (int i = 1; i <= 3; ++i) s.std_verts[i] = StdVertex{120.0 * (i - 1), lat, true};
  // Vertex 1 on the prime meridian.
  s.default_orientation = Orientation{0.0, 90.0, 180.0};
  // Face 3 inverted in the middle of a triangle of edge 2, the three faces
  // around vertex 0 at its corners.
  s.plane = {
    {2, 1.0, 0.0}, {0, 0.5, kSin60}, {4, 1.5, kSin60}, {3, 1.5, kSin60}
  };
  // tan g = 2 sqrt(2)
  s.el_angle = std::atan(2.0 * std::sqrt(2.0));
  s.g_deg = 60.0;
  return s;
}

// ---- Lattice arithmetic for the derivation ----
// A point of the hexagonal lattice as (i, j): x = i - j / 2, y = j sin 60.
struct Lat { int i, j; };
Lat lat_rot60(Lat a, int n) {
  n = ((n % 6) + 6) % 6;
  for (int k = 0; k < n; ++k) a = {a.i - a.j, a.i};
  return a;
}
Lat lat_add(Lat a, Lat b) { return {a.i + b.i, a.j + b.j}; }
Lat lat_sub(Lat a, Lat b) { return {a.i - b.i, a.j - b.j}; }
// Direction of a lattice vector as a multiple of 60 degrees
int lat_dir60(Lat a) {
  const double x = a.i - 0.5 * a.j, y = a.j * kSin60;
  const int k = static_cast<int>(std::lround(std::atan2(y, x) / kPiOver3));
  return ((k % 6) + 6) % 6;
}
double lat_x(Lat a) { return static_cast<double>(a.i) - 0.5 * static_cast<double>(a.j); }
double lat_y(Lat a) { return static_cast<double>(a.j) * kSin60; }

// Face vertices v0, v1, v2 in triangle coordinates: top, lower left, lower
// right.
const Lat kTriVertex[3] = {{1, 1}, {0, 0}, {1, 0}};

// Start of each quad-plane region, in multiples of 60 degrees
// (compute_subtriangle()): region 1 spans 0-60, region 0 60-120, and so on.
const int kRegionStart[6] = {1, 0, 5, 4, 3, 2};

int local_index(const std::array<int, 3>& f, int v) {
  return (f[0] == v) ? 0 : (f[1] == v) ? 1 : (f[2] == v) ? 2 : -1;
}

void fail(const SolidSpec& s, const std::string& what) {
  throw std::logic_error(std::string(s.name) + " topology: " + what);
}

// The faces around vertex v, counter-clockwise seen from outside, from `first`
std::vector<int> face_fan(const SolidSpec& s, int v, int first) {
  std::vector<int> ring{first};
  for (;;) {
    const auto& f = s.faces[ring.back()];
    const int end_v = f[(local_index(f, v) + 2) % 3];
    int next = -1;
    for (int g = 0; g < static_cast<int>(s.faces.size()); ++g) {
      const int k = local_index(s.faces[g], v);
      if (k >= 0 && s.faces[g][(k + 1) % 3] == end_v) next = g;
    }
    if (next < 0) fail(s, "faces around a vertex do not close");
    if (next == first) return ring;
    ring.push_back(next);
  }
}

SolidTopology derive_topology(SolidSpec s) {
  SolidTopology t{};
  const int V = s.n_verts;
  const int F = static_cast<int>(s.faces.size());
  if (V > kMaxVerts || F > kMaxFaces || V != F / 2 + 2) fail(s, "not a solid of triangles");
  t.solid = s.solid;
  t.name = s.name;
  t.n_faces = F;
  t.n_verts = V;
  t.has_quads = !s.diamonds.empty();
  t.snyder = snyder_params(s.el_angle, s.g_deg, F);
  t.vgc = vertex_gc_params(s.el_angle, s.g_deg);
  t.default_orientation = s.default_orientation;
  for (int v = 0; v < V; ++v) t.std_verts[v] = s.std_verts[v];
  for (int f = 0; f < F; ++f) t.plane[f] = s.plane[f];

  for (int v = 0; v < V; ++v) {
    std::vector<int>& nb = t.neighbors[v];
    for (const auto& f : s.faces) {
      const int k = local_index(f, v);
      if (k >= 0) { nb.push_back(f[(k + 1) % 3]); nb.push_back(f[(k + 2) % 3]); }
    }
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    t.valence[v] = static_cast<int>(nb.size());
  }

  for (int q = 0; q < V; ++q) t.fold_axis[q] = -1;
  t.has_folds = false;
  if (t.has_quads) {
    // Counter-clockwise from the quad's vertex, a diamond's primary face
    // runs to the far corner and then to the j corner, its secondary face to
    // the i corner and then to the far corner. A face's vertex list keeps the
    // order the spec gives it, which fixes the vertex its azimuth is read
    // from; the regions below place it in the quad whichever vertex comes
    // first.
    std::vector<int> face_quad(F, -1), face_sub(F, -1);
    for (int q = 1; q <= V - 2; ++q) {
      const int p = s.diamonds[q - 1][0], c = s.diamonds[q - 1][1];
      const int kp = local_index(s.faces[p], q), kc = local_index(s.faces[c], q);
      if (kp < 0 || kc < 0) fail(s, "a diamond face misses its quad's vertex");
      const int far = s.faces[p][(kp + 1) % 3], j_corner = s.faces[p][(kp + 2) % 3];
      const int i_corner = s.faces[c][(kc + 1) % 3];
      if (s.faces[c][(kc + 2) % 3] != far) fail(s, "a diamond's faces do not share its diagonal");
      if (face_quad[p] >= 0 || face_quad[c] >= 0) fail(s, "a face lies in two diamonds");
      face_quad[p] = q; face_sub[p] = 0;
      face_quad[c] = q; face_sub[c] = 1;
      t.corner[q] = {q, i_corner, far, j_corner};
    }
    for (int f = 0; f < F; ++f) if (face_quad[f] < 0) fail(s, "a face lies in no diamond");
    for (int q : {0, V - 1}) t.corner[q] = {q, q, q, q};

    // Regions around each quad's origin
    for (int q = 0; q < V; ++q) {
      for (auto& r : t.region[q]) r = QuadRegion{-1, 0.0, 0.0, 0, false};
      const int first = (q == 0) ? s.pole_first_face[0]
                      : (q == V - 1) ? s.pole_first_face[1]
                      : s.diamonds[q - 1][1];
      const std::vector<int> fan = face_fan(s, q, first);
      const int n = static_cast<int>(fan.size()), n_ccw = s.n_ccw[q];
      if (n > 6 || n_ccw < 1 || n_ccw > n) fail(s, "a quad frame cannot hold its vertex's faces");
      static const int kCcwRegion[6] = {1, 0, 5, 4, 3, 2};
      for (int a = 0; a < n; ++a) {
        const int r = (a < n_ccw) ? kCcwRegion[a] : 2 + (n - 1 - a);
        const int f = fan[a];
        if (t.region[q][r].keep) fail(s, "two faces in one region");
        // The face's vertex counter-clockwise after the origin lies on the
        // region's starting ray.
        const int k = local_index(s.faces[f], q);
        const int beta = lat_dir60(lat_sub(kTriVertex[(k + 1) % 3], kTriVertex[k]));
        const int rot = (((kRegionStart[r] - beta) % 6) + 6) % 6;
        const Lat trans = lat_rot60(kTriVertex[k], rot);
        t.region[q][r] = QuadRegion{f, lat_x(trans), lat_y(trans), rot, true};
      }
    }
    for (int f = 0; f < F; ++f) {
      const QuadRegion& r = t.region[face_quad[f]][face_sub[f] == 0 ? 0 : 1];
      t.placement[f] = FacePlacement{face_quad[f], r.rot60, r.trans_x, r.trans_y};
    }

    // Maps across quad edges. An edge joins two corners of its quad; the face
    // across it lies in another diamond, where the same two vertices are
    // corners too, and the lattice isometry taking one pair to the other
    // carries every point across.
    const Lat corner_pos[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    struct EdgeDef { int a, b; bool cross_i; bool far; };
    const EdgeDef edges[4] = {{kCornerOrigin, kCornerJ, true, false},
                              {kCornerOrigin, kCornerI, false, false},
                              {kCornerI, kCornerFar, true, true},
                              {kCornerJ, kCornerFar, false, true}};
    std::vector<int> fold_edges(V, 0);
    for (int q = 1; q <= V - 2; ++q) {
      for (int e = 0; e < 4; ++e) {
        const int A = t.corner[q][edges[e].a], B = t.corner[q][edges[e].b];
        int nq = -1;
        for (int f = 0; f < F; ++f) {
          if (face_quad[f] != q && local_index(s.faces[f], A) >= 0 &&
              local_index(s.faces[f], B) >= 0) nq = face_quad[f];
        }
        if (nq < 0) fail(s, "no face across a quad edge");
        int ia = -1, ib = -1;
        for (int c = 0; c < 4; ++c) {
          if (t.corner[nq][c] == A) ia = c;
          if (t.corner[nq][c] == B) ib = c;
        }
        if (ia < 0 || ib < 0) fail(s, "a quad edge is no edge of its neighbour");
        const Lat pA = corner_pos[edges[e].a], pB = corner_pos[edges[e].b];
        const Lat qA = corner_pos[ia], qB = corner_pos[ib];
        const int rot = ((lat_dir60(lat_sub(qB, qA)) - lat_dir60(lat_sub(pB, pA))) % 6 + 6) % 6;
        const Lat shift = lat_sub(qA, lat_rot60(pA, rot));
        const Lat cross = edges[e].cross_i ? Lat{1, 0} : Lat{0, 1};
        const Lat along = edges[e].cross_i ? Lat{0, 1} : Lat{1, 0};
        const Lat base = edges[e].far ? cross : Lat{0, 0};
        const Lat k0 = lat_add(lat_rot60(base, rot), shift);
        const Lat ka = lat_rot60(along, rot), kd = lat_rot60(cross, rot);
        QuadEdgeMap& m = t.edge[q][e];
        m.quad = nq;
        m.pole = -1;
        if (edges[e].far && t.is_pole(t.corner[q][edges[e].a])) m.pole = t.corner[q][edges[e].a];
        m.k[0][0] = k0.i; m.k[0][1] = ka.i; m.k[0][2] = kd.i;
        m.k[1][0] = k0.j; m.k[1][1] = ka.j; m.k[1][2] = kd.j;

        // An edge that meets an edge of its own kind: of a near pair the
        // higher quad gives up the interior, of a far pair it takes it.
        const bool nb_near = (ia == kCornerOrigin || ib == kCornerOrigin);
        if (nb_near == !edges[e].far && nq == q) fail(s, "a quad edge folds onto its own quad");
        if (nb_near == !edges[e].far && nq < q) {
          const int axis = edges[e].cross_i ? 0 : 1;
          if (t.fold_axis[q] >= 0 && t.fold_axis[q] != axis) fail(s, "a quad folds along both axes");
          t.fold_axis[q] = axis;
          t.has_folds = true;
          fold_edges[q] += edges[e].far ? 1 : 2;
        }
      }
    }
    for (int q = 1; q <= V - 2; ++q) {
      if (fold_edges[q] != 0 && fold_edges[q] != 3) {
        fail(s, "a quad gives up a near edge without taking the far edge across its box");
      }
    }
  }
  for (int f = 0; f < F; ++f) t.faces[f] = s.faces[f];
  return t;
}

} // anon

const SolidTopology& solid_topology(Solid s) {
  static const SolidTopology icosa = derive_topology(icosahedron_spec());
  static const SolidTopology octa = derive_topology(octahedron_spec());
  static const SolidTopology tetra = derive_topology(tetrahedron_spec());
  switch (s) {
    case Solid::Icosahedron: return icosa;
    case Solid::Octahedron: return octa;
    case Solid::Tetrahedron: return tetra;
  }
  throw std::invalid_argument("unknown solid");
}

// ============================================================================
// Solids on the sphere
// ============================================================================

namespace {

// One table per (solid, orientation) read so far. A few grids are live at
// once, so the list stays short; past kMaxTables the oldest table that is
// neither a default nor the active one is dropped.
constexpr std::size_t kMaxTables = 16;
struct Table { Solid solid; Orientation orientation; std::unique_ptr<PolyData> data; };
std::vector<Table> g_tables;
Orientation g_icosa_default;
const PolyData* g_active = nullptr;

// The face table of solid `t` placed by `o`.
void build_table(const SolidTopology& t, const Orientation& o, PolyData& P) {
  P.topo = &t;
  const int V = t.n_verts, F = t.n_faces;
  const Geo S_pt(deg2rad(o.vert0_lon_deg), deg2rad(o.vert0_lat_deg));
  const double S_az = deg2rad(o.azimuth_deg);
  const Geo newnpold(0.0, S_pt.lat);

  std::array<Geo, kMaxVerts> verts;
  verts[0] = S_pt;
  for (int i = 1; i < V; ++i) {
    const StdVertex& sv = t.std_verts[i];
    const double lon = sv.turns ? wrap_lon(-S_az + deg2rad(sv.lon_deg)) : deg2rad(sv.lon_deg);
    verts[i] = coordtrans(newnpold, Geo(lon, deg2rad(sv.lat_deg)), S_pt.lon);
  }

  for (int i = 0; i < F; ++i) {
    Geo tri[3] = { verts[t.faces[i][0]], verts[t.faces[i][1]], verts[t.faces[i][2]] };
    Geo c = sph_tricen(tri);
    P.centers[i]       = c;
    P.center_sinlat[i] = std::sin(c.lat);
    P.center_coslat[i] = std::cos(c.lat);
    P.center_lon[i]    = c.lon;
  }

  for (int i = 0; i < F; ++i) {
    const Geo& c  = P.centers[i];
    const Geo& t0 = verts[t.faces[i][0]];
    const double num = std::cos(t0.lat) * std::sin(t0.lon - c.lon);
    const double den = P.center_coslat[i] * std::sin(t0.lat)
                     - std::sin(c.lat) * std::cos(t0.lat) * std::cos(t0.lon - c.lon);
    P.face_azimuth_offset[i] = std::atan2(num, den);
  }

  for (int v = 0; v < V; ++v) P.verts[v] = verts[v];

  // A face's corners project to its triangle's corners, so the three pairs
  // (triangle coordinates, vertex position) fix the face's affine map. Every
  // face projection takes a vertex to the same corner; the corners are read
  // through the ISEA projection, so the table does not depend on which
  // projection is active when it is built.
  for (int i = 0; i < F; ++i) {
    Vec3 Pk[3];
    double tk[3][2];
    for (int k = 0; k < 3; ++k) {
      Pk[k] = ll2xyz(verts[t.faces[i][k]]);
      const auto xy = project_to_face(verts[t.faces[i][k]], P, i, FaceProjection::ISEA);
      tk[k][0] = xy.first;
      tk[k][1] = xy.second;
    }
    const double a = tk[1][0] - tk[0][0], b = tk[2][0] - tk[0][0];
    const double c = tk[1][1] - tk[0][1], d = tk[2][1] - tk[0][1];
    const double det = a * d - b * c;
    const double inv[2][2] = { {  d / det, -b / det },
                               { -c / det,  a / det } };
    const double E1[3] = { Pk[1].x - Pk[0].x, Pk[1].y - Pk[0].y, Pk[1].z - Pk[0].z };
    const double E2[3] = { Pk[2].x - Pk[0].x, Pk[2].y - Pk[0].y, Pk[2].z - Pk[0].z };
    const double P0[3] = { Pk[0].x, Pk[0].y, Pk[0].z };
    for (int r = 0; r < 3; ++r) {
      P.solid_x[i][r] = E1[r] * inv[0][0] + E2[r] * inv[1][0];
      P.solid_y[i][r] = E1[r] * inv[0][1] + E2[r] * inv[1][1];
      P.solid_origin[i][r] = P0[r] - P.solid_x[i][r] * tk[0][0]
                                   - P.solid_y[i][r] * tk[0][1];
    }

    // Edge planes, each normal pointing towards the face's centre
    const Vec3 ctr{ Pk[0].x + Pk[1].x + Pk[2].x, Pk[0].y + Pk[1].y + Pk[2].y,
                    Pk[0].z + Pk[1].z + Pk[2].z };
    for (int k = 0; k < 3; ++k) {
      const Vec3& u = Pk[k];
      const Vec3& w = Pk[(k + 1) % 3];
      double n[3] = { u.y * w.z - u.z * w.y, u.z * w.x - u.x * w.z, u.x * w.y - u.y * w.x };
      const double inv_len = 1.0 / std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      for (int r = 0; r < 3; ++r) n[r] = inv_len * n[r];
      const bool flip = n[0] * ctr.x + n[1] * ctr.y + n[2] * ctr.z < 0.0;
      for (int r = 0; r < 3; ++r) P.edge_normal[i][k][r] = flip ? -1.0 * n[r] : n[r];
    }
  }
}

} // anon

Orientation default_orientation(Solid s) {
  if (s == Solid::Icosahedron) return g_icosa_default;
  return solid_topology(s).default_orientation;
}

void use_solid(Solid s, const Orientation& o) {
  if (!std::isfinite(o.vert0_lon_deg) || !std::isfinite(o.vert0_lat_deg) ||
      !std::isfinite(o.azimuth_deg)) {
    throw std::invalid_argument("orientation: vert0_lon, vert0_lat and azimuth must be finite");
  }
  if (o.vert0_lat_deg < -90.0 || o.vert0_lat_deg > 90.0) {
    throw std::invalid_argument("orientation: vert0_lat must lie in [-90, 90]");
  }
  for (const Table& t : g_tables) {
    if (t.solid == s && t.orientation == o) {
      g_active = t.data.get();
      return;
    }
  }
  if (g_tables.size() >= kMaxTables) {
    for (auto it = g_tables.begin(); it != g_tables.end(); ++it) {
      const bool is_default = (it->orientation == default_orientation(it->solid));
      if (!is_default && it->data.get() != g_active) {
        g_tables.erase(it);
        break;
      }
    }
  }
  std::unique_ptr<PolyData> data(new PolyData());
  build_table(solid_topology(s), o, *data);
  g_active = data.get();
  g_tables.push_back(Table{s, o, std::move(data)});
}

void use_default_solid(Solid s) { use_solid(s, default_orientation(s)); }

void use_default_orientation() { use_default_solid(Solid::Icosahedron); }

void build_icosa_full(double vert0_lon_deg, double vert0_lat_deg, double azimuth_deg) {
  Orientation o;
  o.vert0_lon_deg = vert0_lon_deg;
  o.vert0_lat_deg = vert0_lat_deg;
  o.azimuth_deg = azimuth_deg;
  use_solid(Solid::Icosahedron, o);
  g_icosa_default = o;
}

const PolyData& poly() {
  if (g_active == nullptr) use_default_orientation();
  return *g_active;
}

const std::array<Geo, kMaxFaces>& face_centers() { return poly().centers; }

double get_face_azimuth_offset(int face) {
  const PolyData& P = poly();
  if (face < 0 || face >= P.n_faces()) return 0.0;
  return P.face_azimuth_offset[face];
}

void face_tri_to_solid(int face, double tx, double ty, double out[3]) {
  const PolyData& S = poly();
  for (int r = 0; r < 3; ++r) {
    out[r] = S.solid_origin[face][r] + tx * S.solid_x[face][r] + ty * S.solid_y[face][r];
  }
}

void face_tri_to_sphere(int face, double tx, double ty, double out[3]) {
  const auto ll = face_xy_to_sphere_ll(tx, ty, face);
  const Vec3 v = ll2xyz(Geo(deg2rad(ll.first), deg2rad(ll.second)));
  out[0] = v.x;
  out[1] = v.y;
  out[2] = v.z;
}

void face_tri_to_plane(int face, double tx, double ty, double& px, double& py) {
  const PlaneTriLayout& layout = topo().plane[face];
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
  const PolyData& P = poly();
  const Geo point(deg2rad(lon_deg), to_sphere_lat(deg2rad(lat_deg)));
  int best = 0;
  double bestd = std::acos(clampd(std::sin(P.centers[0].lat)*std::sin(point.lat)
                       + std::cos(P.centers[0].lat)*std::cos(point.lat)*std::cos(P.centers[0].lon - point.lon),
                       -1.0, 1.0));
  for (int i = 1; i < P.n_faces(); ++i) {
    const auto& c = P.centers[i];
    const double cc = clampd(std::sin(c.lat)*std::sin(point.lat)
                     + std::cos(c.lat)*std::cos(point.lat)*std::cos(c.lon - point.lon), -1.0, 1.0);
    const double d = std::acos(cc);
    if (d < bestd) { best = i; bestd = d; }
  }

  return best;
}

} // namespace hexify
