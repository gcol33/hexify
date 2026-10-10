// hex9_warp.cpp - the trained warp of Hex9's own projection (see hex9_warp.h)
//
// Written from Griffin (2026), Sec. 11b, the .h9warp file format and the
// Clough-Tocher evaluation of libhex9 (core/h9_warp_io.h, h9_warp_fund.h,
// h9_ct.h, h9_warp.h; https://github.com/MrBenGriffin/libhex9), used as a
// reference. The field values are libhex9's (Copyright (c) 2025-2026 Ben
// Griffin, Apache License 2.0), read at run time from the user's copy.
//
// The field lives in libhex9's chart of a face: a triangle of edge sqrt(2)
// centred on the origin, pointing down, its lower corner C at (0, -2H/3) and
// its top edge at y = H/3, H = sqrt(6)/2. Its level-L lattice (n = 3^(L+1)
// steps along an edge) has its points at (i u_x, -2H/3 + j u_y),
// u_x = sqrt(2) / (2n), u_y = H / n, |i| <= j <= n, i + j even. The file
// holds the points of the wedge 0 <= i <= j, i + 3 j <= 2n -- the sixth of
// the face between the line x = 0, the median from the centre to the middle
// of the right edge, and that edge -- in the order of increasing i, then j,
// six numbers each: the displacement (dx, dy) and its gradients
// (d dx/dx, d dx/dy, d dy/dx, d dy/dy). The field elsewhere follows from its
// symmetry: d(T p) = T d(p) for the six symmetries T of the face, and across
// the face's edge by the reflection there.

#include "hex9_warp.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace hexify {
namespace hex9 {

namespace {

struct Field {
  bool ready = false;
  int level = 0;
  int64_t n = 0;                  // lattice steps along an edge
  double ux = 0, uy = 0;          // lattice pitch
  std::vector<double> rows;       // six numbers per wedge point
  std::vector<int64_t> col_start; // first row of each wedge column i
  std::vector<int64_t> col_j0;    // its first j
};

Field g_field;

const double kW = std::sqrt(2.0);
const double kH = kW * std::sqrt(3.0) * 0.5;
const double kVF = -(kH / 3.0) * 2.0;     // lower corner C
const double kR3 = std::sqrt(3.0);

// A 2 x 2 matrix, row-major
struct Mat { double a, b, c, d; };

inline Mat mul(const Mat& p, const Mat& q) {
  return {p.a * q.a + p.b * q.c, p.a * q.b + p.b * q.d,
          p.c * q.a + p.d * q.c, p.c * q.b + p.d * q.d};
}
inline Mat transpose(const Mat& p) { return {p.a, p.c, p.b, p.d}; }

// The face's symmetries in libhex9's fold order (identity, the two turns,
// then each of them after the mirror x -> -x), the wedge's median normal and
// the right edge's reflection about C.
struct Geometry {
  Mat T[6];
  double n2x, n2y;       // wedge side of the median: p . n2 >= 0
  double n3x, n3y;       // face side of the right edge: (p - C) . n3 >= 0
  Mat s3;                // reflection across the right edge
};

Geometry make_geometry() {
  Geometry g;
  const double Cx = 0.0, Cy = kVF, Bx = kW / 2.0, By = kH / 3.0;
  const double Mx = 0.5 * (Cx + Bx), My = 0.5 * (Cy + By);
  const double mn = std::sqrt(Mx * Mx + My * My);
  const double u2x = Mx / mn, u2y = My / mn;
  g.n2x = -u2y;
  g.n2y = u2x;
  if (g.n2x * Cx + g.n2y * Cy < 0.0) { g.n2x = -g.n2x; g.n2y = -g.n2y; }
  const double ex = Bx - Cx, ey = By - Cy;
  const double en = std::sqrt(ex * ex + ey * ey);
  const double u3x = ex / en, u3y = ey / en;
  g.n3x = -u3y;
  g.n3y = u3x;
  if (g.n3x * (0.0 - Cx) + g.n3y * (0.0 - Cy) < 0.0) { g.n3x = -g.n3x; g.n3y = -g.n3y; }
  g.s3 = {2.0 * u3x * u3x - 1.0, 2.0 * u3x * u3y, 2.0 * u3x * u3y, 2.0 * u3y * u3y - 1.0};
  const Mat base[3] = {{1.0, 0.0, 0.0, 1.0},
                       {-0.5, -kR3 / 2.0, kR3 / 2.0, -0.5},
                       {-0.5, kR3 / 2.0, -kR3 / 2.0, -0.5}};
  for (int i = 0; i < 3; i++) {
    g.T[i] = base[i];
    g.T[3 + i] = {-base[i].a, -base[i].b, base[i].c, base[i].d};
  }
  return g;
}

const Geometry& geometry() {
  static const Geometry g = make_geometry();
  return g;
}

// The symmetry carrying (x, y) into the wedge's cone, and the image
int fold(const Geometry& g, double x, double y, double& fx, double& fy) {
  for (int i = 0; i < 6; i++) {
    const Mat& T = g.T[i];
    const double cx = x * T.a + y * T.b, cy = x * T.c + y * T.d;
    if (cx >= -1e-9 && cx * g.n2x + cy * g.n2y >= -1e-9) {
      fx = cx;
      fy = cy;
      return i;
    }
  }
  fx = x;
  fy = y;
  return 0;
}

// A field value and its gradients at one lattice point
struct Datum {
  double d[2];     // dx, dy
  double J[2][2];  // J[r][c] = d d_r / d x_c
};

// The datum of the wedge point nearest (x, y), which must be one
bool wedge_datum(const Field& f, double x, double y, Datum& out) {
  const int64_t i = std::llround(x / f.ux);
  const int64_t j = std::llround((y - kVF) / f.uy);
  if (i < 0 || i > f.n || j < f.col_j0[i] || ((i + j) & 1) || i > j ||
      i + 3 * j > 2 * f.n) {
    return false;
  }
  const double* r = &f.rows[6 * static_cast<size_t>(f.col_start[i] + (j - f.col_j0[i]) / 2)];
  out.d[0] = r[0];
  out.d[1] = r[1];
  out.J[0][0] = r[2];
  out.J[0][1] = r[3];
  out.J[1][0] = r[4];
  out.J[1][1] = r[5];
  return true;
}

// The datum at lattice point (x, y) of the face or its edge band: carried
// into the wedge by the face's symmetries and, past the right edge, the
// reflection there. With the chain p = M q, d(p) = M d_w(q) and
// J(p) = M J_w M^T.
Datum lattice_datum(const Field& f, double x, double y) {
  const Geometry& g = geometry();
  Mat M = {1.0, 0.0, 0.0, 1.0};
  for (int pass = 0; pass < 3; pass++) {
    double qx, qy;
    const int k = fold(g, x, y, qx, qy);
    M = mul(M, transpose(g.T[k]));
    if ((qx - 0.0) * g.n3x + (qy - kVF) * g.n3y >= -1e-9) {
      Datum w;
      if (wedge_datum(f, qx, qy, w)) {
        Datum out;
        const Mat J = {w.J[0][0], w.J[0][1], w.J[1][0], w.J[1][1]};
        const Mat MJ = mul(mul(M, J), transpose(M));
        out.d[0] = M.a * w.d[0] + M.b * w.d[1];
        out.d[1] = M.c * w.d[0] + M.d * w.d[1];
        out.J[0][0] = MJ.a;
        out.J[0][1] = MJ.b;
        out.J[1][0] = MJ.c;
        out.J[1][1] = MJ.d;
        return out;
      }
      break;
    }
    // Past the right edge: reflect across it
    const double rx = qx - 0.0, ry = qy - kVF;
    x = 0.0 + g.s3.a * rx + g.s3.b * ry;
    y = kVF + g.s3.c * rx + g.s3.d * ry;
    M = mul(M, g.s3);
  }
  Datum zero = {{0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}}};
  return zero;
}

// The Clough-Tocher cubic of one field component over a triangle P[3], from
// its values v and gradients gr at the corners, at barycentric b: the
// triangle is split at its centroid into three cubics joined C1 (Alfeld),
// as scipy's CloughTocher2DInterpolator builds them, the cross-boundary
// derivatives taken from the centroid of the triangle across each edge.
double clough_tocher(const double P[3][2], const double v[3], const double gr[3][2],
                     const double G[3], const double b[3]) {
  const double e01x = P[1][0] - P[0][0], e01y = P[1][1] - P[0][1];
  const double e02x = P[2][0] - P[0][0], e02y = P[2][1] - P[0][1];
  const double e12x = P[2][0] - P[1][0], e12y = P[2][1] - P[1][1];
  const double d01 = gr[0][0] * e01x + gr[0][1] * e01y;
  const double d02 = gr[0][0] * e02x + gr[0][1] * e02y;
  const double d10 = -(gr[1][0] * e01x + gr[1][1] * e01y);
  const double d12 = gr[1][0] * e12x + gr[1][1] * e12y;
  const double d20 = -(gr[2][0] * e02x + gr[2][1] * e02y);
  const double d21 = -(gr[2][0] * e12x + gr[2][1] * e12y);

  const double c3000 = v[0], c0300 = v[1], c0030 = v[2];
  const double c2100 = (d01 + 3.0 * c3000) / 3.0;
  const double c1200 = (d10 + 3.0 * c0300) / 3.0;
  const double c2010 = (d02 + 3.0 * c3000) / 3.0;
  const double c0210 = (d12 + 3.0 * c0300) / 3.0;
  const double c1020 = (d20 + 3.0 * c0030) / 3.0;
  const double c0120 = (d21 + 3.0 * c0030) / 3.0;
  const double c2001 = (c2100 + c2010 + c3000) / 3.0;
  const double c0201 = (c1200 + c0300 + c0210) / 3.0;
  const double c0021 = (c1020 + c0120 + c0030) / 3.0;
  const double c0111 = (G[0] * (-c0300 + 3.0 * c0210 - 3.0 * c0120 + c0030)
                        + (-c0300 + 2.0 * c0210 - c0120 + c0021 + c0201)) / 2.0;
  const double c1011 = (G[1] * (-c0030 + 3.0 * c1020 - 3.0 * c2010 + c3000)
                        + (-c0030 + 2.0 * c1020 - c2010 + c2001 + c0021)) / 2.0;
  const double c1101 = (G[2] * (-c3000 + 3.0 * c2100 - 3.0 * c1200 + c0300)
                        + (-c3000 + 2.0 * c2100 - c1200 + c2001 + c0201)) / 2.0;
  const double c1002 = (c1101 + c1011 + c2001) / 3.0;
  const double c0102 = (c1101 + c0111 + c0201) / 3.0;
  const double c0012 = (c1011 + c0111 + c0021) / 3.0;
  const double c0003 = (c1002 + c0102 + c0012) / 3.0;

  const double mn = std::min(b[0], std::min(b[1], b[2]));
  const double s1 = b[0] - mn, s2 = b[1] - mn, s3 = b[2] - mn, s4 = 3.0 * mn;
  return s1 * s1 * s1 * c3000 + 3.0 * s1 * s1 * s2 * c2100 + 3.0 * s1 * s1 * s3 * c2010 +
         3.0 * s1 * s1 * s4 * c2001 + 3.0 * s1 * s2 * s2 * c1200 +
         6.0 * s1 * s2 * s4 * c1101 + 3.0 * s1 * s3 * s3 * c1020 + 6.0 * s1 * s3 * s4 * c1011 +
         3.0 * s1 * s4 * s4 * c1002 + s2 * s2 * s2 * c0300 + 3.0 * s2 * s2 * s3 * c0210 +
         3.0 * s2 * s2 * s4 * c0201 + 3.0 * s2 * s3 * s3 * c0120 + 6.0 * s2 * s3 * s4 * c0111 +
         3.0 * s2 * s4 * s4 * c0102 + s3 * s3 * s3 * c0030 + 3.0 * s3 * s3 * s4 * c0021 +
         3.0 * s3 * s4 * s4 * c0012 + s4 * s4 * s4 * c0003;
}

// The displacement at chart point (x, y) of the wedge's cone, from the
// lattice triangle holding it
void wedge_delta(const Field& f, double x, double y, double& dx, double& dy) {
  // Lattice coordinates along the 60- and 120-degree steps
  const double s = x / f.ux, t = (y - kVF) / f.uy;
  const double u = 0.5 * (s + t), v = 0.5 * (t - s);
  const double U = std::floor(u), V = std::floor(v);
  // The unit rhombus at (U, V) splits into the triangle on its lower side
  // and the one on its upper side.
  const bool lower = (u - U) + (v - V) < 1.0;
  const double cu[3] = {lower ? U : U + 1.0, lower ? U + 1.0 : U, lower ? U : U + 1.0};
  const double cv[3] = {V, lower ? V : V + 1.0, V + 1.0};
  double P[3][2];
  for (int k = 0; k < 3; k++) {
    P[k][0] = (cu[k] - cv[k]) * f.ux;
    P[k][1] = kVF + (cu[k] + cv[k]) * f.uy;
  }
  Datum D[3];
  for (int k = 0; k < 3; k++) D[k] = lattice_datum(f, P[k][0], P[k][1]);

  // Barycentric coordinates of the point
  const double T00 = P[0][0] - P[2][0], T01 = P[1][0] - P[2][0];
  const double T10 = P[0][1] - P[2][1], T11 = P[1][1] - P[2][1];
  const double det = T00 * T11 - T01 * T10;
  const double ex = x - P[2][0], ey = y - P[2][1];
  double b[3];
  b[0] = (T11 * ex - T01 * ey) / det;
  b[1] = (-T10 * ex + T00 * ey) / det;
  b[2] = 1.0 - b[0] - b[1];

  // Cross-boundary weights from the centroid of the triangle across each
  // edge, the reflection of this one's third corner through the edge
  const double V4x = (P[0][0] + P[1][0] + P[2][0]) / 3.0;
  const double V4y = (P[0][1] + P[1][1] + P[2][1]) / 3.0;
  double G[3];
  for (int k = 0; k < 3; k++) {
    const double* A = P[(k + 1) % 3];
    const double* B = P[(k + 2) % 3];
    const double nx = (2.0 * A[0] + 2.0 * B[0] - P[k][0]) / 3.0;
    const double ny = (2.0 * A[1] + 2.0 * B[1] - P[k][1]) / 3.0;
    const double ddx = nx - V4x, ddy = ny - V4y;
    const double ax = V4x - A[0], ay = V4y - A[1];
    const double bx = B[0] - A[0], by = B[1] - A[1];
    G[k] = (ddy * ax - ddx * ay) / (ddx * by - ddy * bx);
  }

  for (int r = 0; r < 2; r++) {
    const double vals[3] = {D[0].d[r], D[1].d[r], D[2].d[r]};
    const double grs[3][2] = {{D[0].J[r][0], D[0].J[r][1]},
                              {D[1].J[r][0], D[1].J[r][1]},
                              {D[2].J[r][0], D[2].J[r][1]}};
    (r == 0 ? dx : dy) = clough_tocher(P, vals, grs, G, b);
  }
}

// The displacement at chart point (x, y): folded into the wedge's cone,
// evaluated there and unfolded, d(p) = T^T d(T p).
void chart_delta(double x, double y, double& dx, double& dy) {
  const Geometry& g = geometry();
  double fx, fy;
  const int k = fold(g, x, y, fx, fy);
  double wx, wy;
  wedge_delta(g_field, fx, fy, wx, wy);
  const Mat& T = g.T[k];
  dx = wx * T.a + wy * T.c;
  dy = wx * T.b + wy * T.d;
}

// Face triangle coordinates to libhex9's chart and back: a turn of the up
// triangle to the down one, scaled from unit edge to sqrt(2)
const double kInR = 0.28867513459481288225;   // 1 / (2 sqrt(3)), the inradius
inline void face_to_chart(double x, double y, double& cx, double& cy) {
  cx = kW * (x - 0.5);
  cy = -kW * (y - kInR);
}
inline void chart_to_face(double cx, double cy, double& x, double& y) {
  x = cx / kW + 0.5;
  y = -cy / kW + kInR;
}

uint32_t crc32(const unsigned char* data, std::size_t len) {
  static uint32_t table[256];
  static bool init = false;
  if (!init) {
    for (uint32_t i = 0; i < 256; i++) {
      uint32_t c = i;
      for (int k = 0; k < 8; k++) c = (c >> 1) ^ ((c & 1u) ? 0xedb88320u : 0u);
      table[i] = c;
    }
    init = true;
  }
  uint32_t c = 0xffffffffu;
  for (std::size_t i = 0; i < len; i++) c = (c >> 8) ^ table[(c ^ data[i]) & 0xffu];
  return c ^ 0xffffffffu;
}

template <typename T>
T read_le(const unsigned char* p) {
  // Little-endian bytes into T, independent of the host's byte order
  unsigned char bytes[sizeof(T)];
  const uint16_t probe = 1;
  const bool little = *reinterpret_cast<const unsigned char*>(&probe) == 1;
  for (std::size_t k = 0; k < sizeof(T); k++) bytes[k] = little ? p[k] : p[sizeof(T) - 1 - k];
  T v;
  std::memcpy(&v, bytes, sizeof(T));
  return v;
}

} // anon

bool warp_load(const unsigned char* data, std::size_t len, std::string& err) {
  constexpr std::size_t kHeader = 56;
  if (len < kHeader || std::memcmp(data, "H9WP", 4) != 0) {
    err = "not a .h9warp file (no H9WP header)";
    return false;
  }
  const uint16_t version = read_le<uint16_t>(data + 4);
  const int level = data[6];
  const uint32_t count = read_le<uint32_t>(data + 8);
  const int dtype = data[12];
  const int flags = data[13];
  const uint32_t crc = read_le<uint32_t>(data + 48);
  if (version != 4 || dtype != 0 || (flags & 6) != 6 || (flags & 1) != 0) {
    err = "the .h9warp file is not libhex9's version-4 wedge field with gradients";
    return false;
  }
  if (len - kHeader != static_cast<std::size_t>(count) * 48) {
    err = "the .h9warp file's length does not match its header";
    return false;
  }
  if (crc32(data + kHeader, len - kHeader) != crc) {
    err = "the .h9warp file fails its CRC check";
    return false;
  }
  Field f;
  f.level = level;
  f.n = 1;
  for (int k = 0; k <= level; k++) f.n *= 3;
  f.ux = (kW / static_cast<double>(f.n)) / 2.0;
  f.uy = kH / static_cast<double>(f.n);
  // Wedge columns i = 0 .. n: j from i to (2n - i) / 3, of i's parity
  f.col_start.assign(f.n + 1, 0);
  f.col_j0.assign(f.n + 1, 0);
  int64_t total = 0;
  for (int64_t i = 0; i <= f.n; i++) {
    const int64_t j0 = i;
    int64_t j1 = (2 * f.n - i) / 3;
    if (((j1 - j0) & 1) != 0) j1--;
    f.col_start[i] = total;
    f.col_j0[i] = j0;
    if (j1 >= j0) total += (j1 - j0) / 2 + 1;
  }
  if (total != static_cast<int64_t>(count)) {
    err = "the .h9warp file's point count does not match its level's wedge";
    return false;
  }
  f.rows.resize(static_cast<std::size_t>(count) * 6);
  for (std::size_t k = 0; k < f.rows.size(); k++) {
    f.rows[k] = read_le<double>(data + kHeader + 8 * k);
  }
  f.ready = true;
  g_field = std::move(f);
  return true;
}

bool warp_ready() { return g_field.ready; }

std::vector<uint32_t> warp_texture_words() {
  std::vector<uint32_t> out;
  if (!g_field.ready) return out;
  out.reserve(2 + g_field.col_start.size() + g_field.rows.size());
  out.push_back(static_cast<uint32_t>(g_field.n));
  for (int64_t s : g_field.col_start) out.push_back(static_cast<uint32_t>(s));
  for (double v : g_field.rows) {
    const float f = static_cast<float>(v);
    uint32_t w;
    std::memcpy(&w, &f, 4);
    out.push_back(w);
  }
  return out;
}

void warp_delta(double x, double y, double& dx, double& dy) {
  double cx, cy, cdx, cdy;
  face_to_chart(x, y, cx, cy);
  chart_delta(cx, cy, cdx, cdy);
  dx = cdx / kW;
  dy = -cdy / kW;
}

void warp_apply(double x, double y, double& px, double& py) {
  double dx, dy;
  warp_delta(x, y, dx, dy);
  px = x + dx;
  py = y + dy;
}

void warp_solve(double px, double py, double& x, double& y) {
  // Newton in the chart, from the point itself, with a difference Jacobian
  // (step 1e-7) as libhex9 solves it
  double tx, ty;
  face_to_chart(px, py, tx, ty);
  double cx = tx, cy = ty;
  const double h = 1e-7;
  for (int it = 0; it < 25; it++) {
    double dx, dy;
    chart_delta(cx, cy, dx, dy);
    const double ex = cx + dx - tx, ey = cy + dy - ty;
    if (std::fabs(ex) + std::fabs(ey) < 1e-15) break;
    double dxx, dyx, dxy, dyy;
    chart_delta(cx + h, cy, dxx, dyx);
    chart_delta(cx, cy + h, dxy, dyy);
    const double a = 1.0 + (dxx - dx) / h, b = (dxy - dx) / h;
    const double c = (dyx - dy) / h, d = 1.0 + (dyy - dy) / h;
    const double det = a * d - b * c;
    if (std::fabs(det) > 1e-12) {
      cx -= (d * ex - b * ey) / det;
      cy -= (a * ey - c * ex) / det;
    } else {
      cx -= ex;
      cy -= ey;
    }
  }
  chart_to_face(cx, cy, x, y);
}

} // namespace hex9
} // namespace hexify
