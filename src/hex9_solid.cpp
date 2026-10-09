// hex9_solid.cpp - Hex9 cells on hexify's octahedron (see hex9_solid.h)

#include "hex9_solid.h"
#include "coordinate_transforms.h"
#include "constants.h"
#include "polyhedron.h"
#include <cmath>
#include <stdexcept>

namespace hexify {
namespace hex9 {

namespace {

// The axis (0 X, 1 Y, 2 Z) and sign of each vertex of the octahedron
const int kVertexAxis[6] = {2, 0, 1, 0, 1, 2};
const int kVertexSign[6] = {1, 1, 1, -1, -1, -1};

const SolidTopology& octahedron() {
  const SolidTopology& t = topo();
  if (t.solid != Solid::Octahedron) {
    throw std::invalid_argument("Hex9 (aperture 9) is defined on the octahedron only");
  }
  return t;
}

// The face holding the octant whose negative axes are the bits of `mask`
int face_of_mask(int mask) {
  static int table[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  if (table[0] < 0) {
    const SolidTopology& t = octahedron();
    for (int f = 0; f < t.n_faces; f++) {
      int m = 0;
      for (int k = 0; k < 3; k++) {
        const int v = t.faces[f][k];
        if (kVertexSign[v] < 0) m |= 1 << kVertexAxis[v];
      }
      table[m] = f;
    }
  }
  return table[mask];
}

// A face's triangle coordinates as weights on its three vertices: the top,
// lower-left and lower-right corners of the unit triangle.
inline void face_weights(double tx, double ty, double b[3]) {
  b[0] = ty / kSin60;
  b[2] = tx - 0.5 * b[0];
  b[1] = 1.0 - b[0] - b[2];
}

} // anon

OctPoint face_point_cell(int face, double tx, double ty, int level) {
  const SolidTopology& t = octahedron();
  double b[3];
  face_weights(tx, ty, b);
  double sum = 0.0;
  for (int k = 0; k < 3; k++) {
    if (!(b[k] > 0.0)) b[k] = 0.0;
    sum += b[k];
  }
  double w[3];
  int s[3];
  for (int k = 0; k < 3; k++) {
    const int v = t.faces[face][k];
    w[kVertexAxis[v]] = b[k] / sum;
    s[kVertexAxis[v]] = kVertexSign[v];
  }
  return locate(w, s, level);
}

void lattice_face_xy(const OctPoint& p, int level, int& face, double& tx, double& ty) {
  const SolidTopology& t = octahedron();
  int mask = 0;
  for (int a = 0; a < 3; a++) if (p.c[a] < 0) mask |= 1 << a;
  face = face_of_mask(mask);
  const double n = static_cast<double>(level_steps(level));
  double b[3];
  for (int k = 0; k < 3; k++) {
    const int64_t c = p.c[kVertexAxis[t.faces[face][k]]];
    b[k] = static_cast<double>(c < 0 ? -c : c) / n;
  }
  tx = 0.5 * b[0] + b[2];
  ty = kSin60 * b[0];
}

void cell_quad_ij(const OctPoint& centre, int level, int& quad, long long& i, long long& j) {
  int face;
  double tx, ty, qx, qy;
  lattice_face_xy(centre, level, face, tx, ty);
  icosa_tri_to_quad_xy(face, tx, ty, quad, qx, qy);
  const long long dim = quad_edge_dim(3, frame_resolution(level));
  const double v = qy / kSin60;
  const double u = qx + 0.5 * v;
  i = std::llround(u * static_cast<double>(dim));
  j = std::llround(v * static_cast<double>(dim));
  if (!substrate_ij_canonicalize(quad, i, j, dim)) {
    throw std::logic_error("hex9: a cell centre owned by no quad");
  }
}

bool quad_ij_lattice(int quad, long long i, long long j, int level, OctPoint& p) {
  const SolidTopology& t = octahedron();
  double qx, qy, tx, ty;
  quad_ij_to_xy(quad, i, j, 3, frame_resolution(level), qx, qy);
  int face;
  if (!try_quad_xy_to_icosa_tri(quad, qx, qy, face, tx, ty)) return false;
  double b[3];
  face_weights(tx, ty, b);
  const int64_t n = level_steps(level);
  int64_t sum = 0;
  int64_t a[3];
  for (int k = 0; k < 3; k++) {
    const double x = b[k] * static_cast<double>(n);
    a[k] = std::llround(x);
    if (std::fabs(x - static_cast<double>(a[k])) > 1e-3 || a[k] < 0) return false;
    sum += a[k];
  }
  if (sum != n) return false;
  for (int k = 0; k < 3; k++) {
    const int v = t.faces[face][k];
    p.c[kVertexAxis[v]] = kVertexSign[v] * a[k];
  }
  return true;
}

} // namespace hex9
} // namespace hexify
