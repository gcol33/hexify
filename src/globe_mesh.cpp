// globe_mesh.cpp
// Surfaces, land and lines placed on the solid and on the sphere
//
// Every point here carries two positions: on the flat face of the solid it
// lies on, and on the unit sphere, through the face projection of that face.
// A renderer can then fold the solid into the sphere by blending the two.
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "globe_mesh.h"
#include "constants.h"
#include "coordinate_transforms.h"
#include "polyhedron.h"
#include "projection_forward.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <deque>
#include <limits>
#include <string>
#include "rcpp_icosa.h"

using namespace Rcpp;

namespace hexify {

// ============================================================================
// FaceMesh
// ============================================================================

size_t FaceMesh::KeyHash::operator()(const Key& k) const {
  uint64_t a, b;
  std::memcpy(&a, &k.tx, sizeof a);
  std::memcpy(&b, &k.ty, sizeof b);
  uint64_t h = a * 0x9E3779B97F4A7C15ULL;
  h ^= b + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  h ^= (static_cast<uint64_t>(k.item) << 5) ^ static_cast<uint64_t>(k.face);
  return static_cast<size_t>(h);
}

int FaceMesh::vertex(int item, int face, double tx, double ty,
                     const double* sphere) {
  Key key{item, face, tx, ty};
  auto it = index_.find(key);
  if (it != index_.end()) return it->second;
  int id = static_cast<int>(item_.size());
  index_.emplace(key, id);
  item_.push_back(item);
  face_.push_back(face);
  tx_.push_back(tx);
  ty_.push_back(ty);
  for (int r = 0; r < 3; r++) {
    sphere_.push_back(sphere ? sphere[r] : NA_REAL);
  }
  return id;
}

void FaceMesh::triangle(int a, int b, int c) {
  if (a == b || b == c || c == a) return;
  tri_.insert(tri_.end(), {a, b, c});
}

void FaceMesh::convex_polygon(int item, int face, const std::vector<double>& tx,
                              const std::vector<double>& ty) {
  size_t n = tx.size();
  if (n < 3) return;
  int first = vertex(item, face, tx[0], ty[0]);
  int prev = vertex(item, face, tx[1], ty[1]);
  for (size_t k = 2; k < n; k++) {
    int cur = vertex(item, face, tx[k], ty[k]);
    triangle(first, prev, cur);
    prev = cur;
  }
}

void FaceMesh::refine(double max_len) {
  if (!(max_len > 0.0)) Rcpp::stop("max_len must be positive");
  double max2 = max_len * max_len;
  // An edge is halved at most once per round and a face edge is about 1
  // long, so the rounds needed grow with log2(1 / max_len).
  for (int round = 0; round < 64; round++) {
    std::unordered_map<uint64_t, int> mid;
    auto midpoint = [&](int a, int b) {
      double dx = tx_[a] - tx_[b], dy = ty_[a] - ty_[b];
      if (dx * dx + dy * dy <= max2) return -1;
      uint64_t key = (static_cast<uint64_t>(std::min(a, b)) << 32) |
                     static_cast<uint32_t>(std::max(a, b));
      auto it = mid.find(key);
      if (it != mid.end()) return it->second;
      int m;
      if (great_circle_) {
        double s[3];
        double len = 0.0;
        for (int r = 0; r < 3; r++) {
          s[r] = sphere_[3 * a + r] + sphere_[3 * b + r];
          len += s[r] * s[r];
        }
        len = std::sqrt(len);
        for (int r = 0; r < 3; r++) s[r] /= len;
        Geo g(std::atan2(s[1], s[0]), std::atan2(s[2], std::hypot(s[0], s[1])));
        auto t = project_to_face(g, poly(), face_[a]);
        m = vertex(item_[a], face_[a], t.first, t.second, s);
      } else {
        m = vertex(item_[a], face_[a], 0.5 * (tx_[a] + tx_[b]),
                   0.5 * (ty_[a] + ty_[b]));
      }
      mid.emplace(key, m);
      return m;
    };

    std::vector<int> out;
    out.reserve(tri_.size());
    bool split = false;
    for (size_t t = 0; t < tri_.size(); t += 3) {
      int v[3] = {tri_[t], tri_[t + 1], tri_[t + 2]};
      int m[3];   // m[k] splits the edge v[k] -> v[k + 1]
      int n_split = 0;
      for (int k = 0; k < 3; k++) {
        m[k] = midpoint(v[k], v[(k + 1) % 3]);
        if (m[k] >= 0) n_split++;
      }
      if (n_split == 0) {
        out.insert(out.end(), {v[0], v[1], v[2]});
        continue;
      }
      split = true;
      if (n_split == 3) {
        out.insert(out.end(), {v[0], m[0], m[2], m[0], v[1], m[1],
                               m[2], m[1], v[2], m[0], m[1], m[2]});
        continue;
      }
      // Turn the triangle so that its first edge is split and, with two
      // splits, its second as well.
      int r = 0;
      if (n_split == 1) {
        while (m[r] < 0) r++;
      } else {
        while (!(m[r] >= 0 && m[(r + 1) % 3] >= 0)) r++;
      }
      int a = v[r], b = v[(r + 1) % 3], c = v[(r + 2) % 3];
      int mab = m[r], mbc = m[(r + 1) % 3];
      if (n_split == 1) {
        out.insert(out.end(), {a, mab, c, mab, b, c});
      } else {
        out.insert(out.end(), {a, mab, c, mab, b, mbc, mab, mbc, c});
      }
    }
    tri_.swap(out);
    if (!split) return;
  }
  Rcpp::stop("mesh refinement did not converge");
}

List FaceMesh::to_list() const {
  size_t n = item_.size();
  NumericVector solid(3 * n), sphere(3 * n);
  IntegerVector item(n);
  double p[3];
  for (size_t k = 0; k < n; k++) {
    face_tri_to_solid(face_[k], tx_[k], ty_[k], p);
    for (int r = 0; r < 3; r++) solid[3 * k + r] = p[r];
    if (ISNA(sphere_[3 * k])) {
      face_tri_to_sphere(face_[k], tx_[k], ty_[k], p);
    } else {
      for (int r = 0; r < 3; r++) p[r] = sphere_[3 * k + r];
    }
    for (int r = 0; r < 3; r++) sphere[3 * k + r] = p[r];
    item[k] = item_[k];
  }
  NumericVector tri(2 * n);
  for (size_t k = 0; k < n; k++) {
    tri[2 * k] = tx_[k];
    tri[2 * k + 1] = ty_[k];
  }
  IntegerVector index(tri_.begin(), tri_.end());
  return List::create(_["solid"] = solid, _["sphere"] = sphere,
                      _["index"] = index, _["item"] = item, _["tri"] = tri);
}

void face_tri_corners(int face, double tx[3], double ty[3]) {
  const PolyData& S = poly();
  for (int k = 0; k < 3; k++) {
    const auto t = project_to_face(S.verts[S.topo->faces[face][k]], S, face);
    tx[k] = t.first;
    ty[k] = t.second;
  }
}

// ============================================================================
// Faces as regions of the sphere
// ============================================================================

namespace {

struct V3 {
  double x, y, z;
};

inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator*(double s, V3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(V3 a, V3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline V3 unit(V3 a) { return (1.0 / std::sqrt(dot(a, a))) * a; }

inline V3 lonlat_vec(double lon_deg, double lat_deg) {
  double lon = lon_deg * kDegToRad, lat = lat_deg * kDegToRad;
  return {std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon),
          std::sin(lat)};
}

inline Geo vec_geo(V3 v) {
  return Geo(std::atan2(v.y, v.x),
             std::atan2(v.z, std::sqrt(v.x * v.x + v.y * v.y)));
}

// Each face of the spherical solid is the region on the inner side of the
// three great circles through its edges: n . p >= 0 for the three normals n
// (PolyData::edge_normal), each pointing towards the face's centre.
inline std::array<V3, 3> face_normals(int f) {
  const auto& e = poly().edge_normal[f];
  return {V3{e[0][0], e[0][1], e[0][2]}, V3{e[1][0], e[1][1], e[1][2]},
          V3{e[2][0], e[2][1], e[2][2]}};
}

// How far inside face f the point p lies: the least of its three plane
// distances, negative outside.
inline double face_margin(int f, V3 p) {
  const auto n = face_normals(f);
  return std::min(dot(n[0], p), std::min(dot(n[1], p), dot(n[2], p)));
}

// The face a point of the sphere lies on; a point on an edge goes to the face
// it lies deepest in, which is either.
int sphere_face(V3 p) {
  int best = 0;
  double best_m = -std::numeric_limits<double>::infinity();
  for (int f = 0; f < poly().n_faces(); f++) {
    double m = face_margin(f, p);
    if (m > best_m) {
      best_m = m;
      best = f;
    }
  }
  return best;
}

inline std::pair<double, double> face_tri(int face, V3 p) {
  return project_to_face(vec_geo(p), poly(), face);
}

// ============================================================================
// Ear slicing of a polygon with holes
// ============================================================================
//
// Triangulates a simple polygon with holes by cutting ears from a linked ring
// of its vertices, after joining each hole to the outer ring by a bridge.
// Rings that touch themselves or cross slightly, as simplified coastlines
// can, are handled the way Mapbox's earcut handles them: by removing repeated
// and collinear points, cutting small self-crossings, and as a last resort
// splitting the ring along a valid diagonal.

class EarSlicer {
public:
  // 'x', 'y': the vertices of all rings, outer ring first; 'ring_start':
  // where each ring starts in them. Appends vertex index triples to 'tri'.
  void run(const std::vector<double>& x, const std::vector<double>& y,
           const std::vector<int>& ring_start, std::vector<int>& tri) {
    tri_ = &tri;
    nodes_.clear();
    int n_ring = static_cast<int>(ring_start.size());
    int total = static_cast<int>(x.size());
    auto ring_end = [&](int r) {
      return r + 1 < n_ring ? ring_start[r + 1] : total;
    };
    Node* outer = linked_list(x, y, ring_start[0], ring_end(0), true);
    if (!outer || outer->next == outer->prev) return;
    if (n_ring > 1) {
      std::vector<Node*> queue;
      for (int r = 1; r < n_ring; r++) {
        Node* list = linked_list(x, y, ring_start[r], ring_end(r), false);
        if (!list) continue;
        if (list == list->next) list->steiner = true;
        queue.push_back(leftmost(list));
      }
      std::sort(queue.begin(), queue.end(), [](Node* a, Node* b) {
        return a->x < b->x || (a->x == b->x && a->y < b->y);
      });
      for (Node* hole : queue) outer = eliminate_hole(hole, outer);
    }
    earcut_linked(outer, 0);
  }

private:
  struct Node {
    int i;
    double x, y;
    Node* prev = nullptr;
    Node* next = nullptr;
    bool steiner = false;
  };

  std::deque<Node> nodes_;
  std::vector<int>* tri_ = nullptr;

  Node* new_node(int i, double x, double y) {
    nodes_.push_back(Node{i, x, y});
    return &nodes_.back();
  }

  static double area(const Node* p, const Node* q, const Node* r) {
    return (q->y - p->y) * (r->x - q->x) - (q->x - p->x) * (r->y - q->y);
  }
  static bool equals(const Node* a, const Node* b) {
    return a->x == b->x && a->y == b->y;
  }
  static int sign(double v) { return (v > 0) - (v < 0); }

  Node* insert_node(int i, double x, double y, Node* last) {
    Node* p = new_node(i, x, y);
    if (!last) {
      p->prev = p;
      p->next = p;
    } else {
      p->next = last->next;
      p->prev = last;
      last->next->prev = p;
      last->next = p;
    }
    return p;
  }
  static void remove_node(Node* p) {
    p->next->prev = p->prev;
    p->prev->next = p->next;
  }

  Node* linked_list(const std::vector<double>& x, const std::vector<double>& y,
                    int start, int end, bool outer) {
    double sum = 0.0;
    for (int i = start, j = end - 1; i < end; j = i++) {
      sum += (x[j] - x[i]) * (y[i] + y[j]);
    }
    Node* last = nullptr;
    if (outer == (sum > 0)) {
      for (int i = start; i < end; i++) last = insert_node(i, x[i], y[i], last);
    } else {
      for (int i = end - 1; i >= start; i--) last = insert_node(i, x[i], y[i], last);
    }
    if (last && equals(last, last->next)) {
      remove_node(last);
      last = last->next;
    }
    return last;
  }

  Node* filter_points(Node* start, Node* end = nullptr) {
    if (!start) return start;
    if (!end) end = start;
    Node* p = start;
    bool again;
    do {
      again = false;
      if (!p->steiner && (equals(p, p->next) || area(p->prev, p, p->next) == 0)) {
        remove_node(p);
        p = end = p->prev;
        if (p == p->next) break;
        again = true;
      } else {
        p = p->next;
      }
    } while (again || p != end);
    return end;
  }

  void emit(const Node* a, const Node* b, const Node* c) {
    tri_->insert(tri_->end(), {a->i, b->i, c->i});
  }

  void earcut_linked(Node* ear, int pass) {
    if (!ear) return;
    Node* stop = ear;
    while (ear->prev != ear->next) {
      Node* prev = ear->prev;
      Node* next = ear->next;
      if (is_ear(ear)) {
        emit(prev, ear, next);
        remove_node(ear);
        ear = next->next;
        stop = next->next;
        continue;
      }
      ear = next;
      if (ear == stop) {
        if (pass == 0) {
          earcut_linked(filter_points(ear), 1);
        } else if (pass == 1) {
          ear = cure_local_intersections(filter_points(ear));
          earcut_linked(ear, 2);
        } else {
          split_earcut(ear);
        }
        break;
      }
    }
  }

  static bool point_in_triangle(double ax, double ay, double bx, double by,
                                double cx, double cy, double px, double py) {
    return (cx - px) * (ay - py) >= (ax - px) * (cy - py) &&
           (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
           (bx - px) * (cy - py) >= (cx - px) * (by - py);
  }

  bool is_ear(const Node* ear) const {
    const Node* a = ear->prev;
    const Node* b = ear;
    const Node* c = ear->next;
    if (area(a, b, c) >= 0) return false;
    double x0 = std::min({a->x, b->x, c->x}), x1 = std::max({a->x, b->x, c->x});
    double y0 = std::min({a->y, b->y, c->y}), y1 = std::max({a->y, b->y, c->y});
    const Node* p = c->next;
    while (p != a) {
      if (p->x >= x0 && p->x <= x1 && p->y >= y0 && p->y <= y1 &&
          !(p->x == a->x && p->y == a->y) &&
          point_in_triangle(a->x, a->y, b->x, b->y, c->x, c->y, p->x, p->y) &&
          area(p->prev, p, p->next) >= 0) {
        return false;
      }
      p = p->next;
    }
    return true;
  }

  static bool on_segment(const Node* p, const Node* q, const Node* r) {
    return q->x <= std::max(p->x, r->x) && q->x >= std::min(p->x, r->x) &&
           q->y <= std::max(p->y, r->y) && q->y >= std::min(p->y, r->y);
  }

  static bool intersects(const Node* p1, const Node* q1, const Node* p2,
                         const Node* q2) {
    int o1 = sign(area(p1, q1, p2)), o2 = sign(area(p1, q1, q2));
    int o3 = sign(area(p2, q2, p1)), o4 = sign(area(p2, q2, q1));
    if (o1 != o2 && o3 != o4) return true;
    if (o1 == 0 && on_segment(p1, p2, q1)) return true;
    if (o2 == 0 && on_segment(p1, q2, q1)) return true;
    if (o3 == 0 && on_segment(p2, p1, q2)) return true;
    if (o4 == 0 && on_segment(p2, q1, q2)) return true;
    return false;
  }

  static bool locally_inside(const Node* a, const Node* b) {
    return area(a->prev, a, a->next) < 0
               ? area(a, b, a->next) >= 0 && area(a, a->prev, b) >= 0
               : area(a, b, a->prev) < 0 || area(a, a->next, b) < 0;
  }

  Node* cure_local_intersections(Node* start) {
    Node* p = start;
    do {
      Node* a = p->prev;
      Node* b = p->next->next;
      if (!equals(a, b) && intersects(a, p, p->next, b) &&
          locally_inside(a, b) && locally_inside(b, a)) {
        emit(a, p, b);
        remove_node(p);
        remove_node(p->next);
        p = start = b;
      }
      p = p->next;
    } while (p != start);
    return filter_points(p);
  }

  static bool intersects_polygon(const Node* a, const Node* b) {
    const Node* p = a;
    do {
      if (p->i != a->i && p->next->i != a->i && p->i != b->i &&
          p->next->i != b->i && intersects(p, p->next, a, b)) {
        return true;
      }
      p = p->next;
    } while (p != a);
    return false;
  }

  static bool middle_inside(const Node* a, const Node* b) {
    const Node* p = a;
    bool inside = false;
    double px = 0.5 * (a->x + b->x), py = 0.5 * (a->y + b->y);
    do {
      if (((p->y > py) != (p->next->y > py)) && p->next->y != p->y &&
          (px < (p->next->x - p->x) * (py - p->y) / (p->next->y - p->y) + p->x)) {
        inside = !inside;
      }
      p = p->next;
    } while (p != a);
    return inside;
  }

  static bool is_valid_diagonal(const Node* a, const Node* b) {
    return a->next->i != b->i && a->prev->i != b->i && !intersects_polygon(a, b) &&
           ((locally_inside(a, b) && locally_inside(b, a) && middle_inside(a, b) &&
             (area(a->prev, a, b->prev) != 0 || area(a, b->prev, b) != 0)) ||
            (equals(a, b) && area(a->prev, a, a->next) > 0 &&
             area(b->prev, b, b->next) > 0));
  }

  // Link a to b by two edges, one each way: the ring splits in two, and the
  // node returned starts the second.
  Node* split_polygon(Node* a, Node* b) {
    Node* a2 = new_node(a->i, a->x, a->y);
    Node* b2 = new_node(b->i, b->x, b->y);
    Node* an = a->next;
    Node* bp = b->prev;
    a->next = b;
    b->prev = a;
    a2->next = an;
    an->prev = a2;
    b2->next = a2;
    a2->prev = b2;
    bp->next = b2;
    b2->prev = bp;
    return b2;
  }

  void split_earcut(Node* start) {
    Node* a = start;
    do {
      Node* b = a->next->next;
      while (b != a->prev) {
        if (a->i != b->i && is_valid_diagonal(a, b)) {
          Node* c = split_polygon(a, b);
          a = filter_points(a, a->next);
          c = filter_points(c, c->next);
          earcut_linked(a, 0);
          earcut_linked(c, 0);
          return;
        }
        b = b->next;
      }
      a = a->next;
    } while (a != start);
  }

  static Node* leftmost(Node* start) {
    Node* p = start;
    Node* best = start;
    do {
      if (p->x < best->x || (p->x == best->x && p->y < best->y)) best = p;
      p = p->next;
    } while (p != start);
    return best;
  }

  static bool sector_contains_sector(const Node* m, const Node* p) {
    return area(m->prev, m, p->prev) < 0 && area(p->next, m, m->next) < 0;
  }

  // The outer ring's vertex a hole's leftmost point can be joined to: cast a
  // ray to the left, take the nearest edge it meets, and among the ring's
  // vertices in the triangle that leaves, the one at the shallowest angle.
  static Node* find_hole_bridge(Node* hole, Node* outer) {
    Node* p = outer;
    double hx = hole->x, hy = hole->y;
    double qx = -std::numeric_limits<double>::infinity();
    Node* m = nullptr;
    do {
      if (hy <= p->y && hy >= p->next->y && p->next->y != p->y) {
        double x = p->x + (hy - p->y) * (p->next->x - p->x) / (p->next->y - p->y);
        if (x <= hx && x > qx) {
          qx = x;
          m = p->x < p->next->x ? p : p->next;
          if (x == hx) return m;
        }
      }
      p = p->next;
    } while (p != outer);
    if (!m) return nullptr;

    Node* stop = m;
    double mx = m->x, my = m->y;
    double tan_min = std::numeric_limits<double>::infinity();
    p = m;
    do {
      if (hx >= p->x && p->x >= mx && hx != p->x &&
          point_in_triangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy,
                            p->x, p->y)) {
        double tan = std::fabs(hy - p->y) / (hx - p->x);
        if (locally_inside(p, hole) &&
            (tan < tan_min ||
             (tan == tan_min &&
              (p->x > m->x || (p->x == m->x && sector_contains_sector(m, p)))))) {
          m = p;
          tan_min = tan;
        }
      }
      p = p->next;
    } while (p != stop);
    return m;
  }

  Node* eliminate_hole(Node* hole, Node* outer) {
    Node* bridge = find_hole_bridge(hole, outer);
    if (!bridge) return outer;
    Node* bridge_reverse = split_polygon(bridge, hole);
    filter_points(bridge_reverse, bridge_reverse->next);
    return filter_points(bridge, bridge->next);
  }
};

// The part of a spherical polygon on the inner side of great circles, given
// by their normals: Sutherland-Hodgman clipping. The polygon's edges are
// great-circle arcs shorter than a half circle, so a point of an edge is a
// normalised blend of its ends, and the clip stays exact.
template <size_t N>
void clip_to_planes(const std::array<V3, N>& n, std::vector<V3>& poly) {
  std::vector<V3> out;
  for (size_t k = 0; k < N && !poly.empty(); k++) {
    out.clear();
    size_t m = poly.size();
    for (size_t i = 0; i < m; i++) {
      V3 a = poly[i], b = poly[(i + 1) % m];
      double da = dot(n[k], a), db = dot(n[k], b);
      if (da >= 0.0) out.push_back(a);
      if ((da >= 0.0) != (db >= 0.0)) {
        out.push_back(unit(a + (da / (da - db)) * (b - a)));
      }
    }
    poly.swap(out);
  }
}

// The sectors of a face, each between the face's centre and a piece of its
// boundary: the arcs from the centre to the corners are where Snyder's
// projection creases, and the vertex-oriented projection (IVEA) creases on
// the arcs to the edge midpoints as well, so it has six sectors and the
// others three. Sector k runs from boundary point k to k + 1 of the ring
// corner 0, (midpoint), corner 1, (midpoint), corner 2, (midpoint). A mesh
// cut into sectors has a vertex on every crease it crosses.
int n_face_sectors() {
  return active_projection() == FaceProjection::IVEA ? 6 : 3;
}

// Boundary point k of face f's sector ring (n_face_sectors()), as a unit
// vector
V3 sector_ring_point(int f, int k) {
  const PolyData& S = poly();
  const int n = n_face_sectors();
  auto corner = [&](int i) {
    const Geo& g = S.verts[S.topo->faces[f][i % 3]];
    return lonlat_vec(g.lon * kRadToDeg, g.lat * kRadToDeg);
  };
  if (n == 3) return corner(k);
  if (k % 2 == 0) return corner(k / 2);
  return unit(corner(k / 2) + corner(k / 2 + 1));
}

// Sector k of face f is the region inside three great circles: the face edge
// and the two arcs from the centre.
std::array<V3, 3> sector_normals(int f, int k) {
  const PolyData& S = poly();
  const Geo& gc = S.centers[f];
  const V3 c = lonlat_vec(gc.lon * kRadToDeg, gc.lat * kRadToDeg);
  const V3 v0 = sector_ring_point(f, k);
  const V3 v1 = sector_ring_point(f, (k + 1) % n_face_sectors());
  // Each normal turned to point towards the sector's third corner.
  auto towards = [](V3 n, V3 p) { return dot(n, p) < 0.0 ? (-1.0) * n : n; };
  return {towards(cross(v0, v1), c), towards(cross(c, v0), v1),
          towards(cross(c, v1), v0)};
}

} // namespace

} // namespace hexify

using namespace hexify;

// The faces of the solid as one mesh, refined until no edge is longer than
// 'max_len' (a face edge is 1). Item k is face k - 1. Each face is laid down
// as its sectors (n_face_sectors()), so the mesh has edges along the
// projection's creases.
// [[Rcpp::export]]
List cpp_globe_faces(NumericVector icosa, double max_len) {
  activate_icosa(icosa);
  FaceMesh mesh;
  double tx[3], ty[3];
  for (int f = 0; f < poly().n_faces(); f++) {
    face_tri_corners(f, tx, ty);
    const double cx = (tx[0] + tx[1] + tx[2]) / 3.0;
    const double cy = (ty[0] + ty[1] + ty[2]) / 3.0;
    // The sector ring in triangle coordinates: corners, with the edge
    // midpoints between them on six sectors
    const int n = n_face_sectors();
    double rx[6], ry[6];
    for (int k = 0; k < n; k++) {
      const int i = n == 3 ? k : k / 2;
      const int i1 = (i + 1) % 3;
      const bool mid = n == 6 && k % 2 == 1;
      rx[k] = mid ? 0.5 * (tx[i] + tx[i1]) : tx[i];
      ry[k] = mid ? 0.5 * (ty[i] + ty[i1]) : ty[i];
    }
    for (int k = 0; k < n; k++) {
      const int k1 = (k + 1) % n;
      mesh.convex_polygon(f + 1, f, {cx, rx[k], rx[k1]}, {cy, ry[k], ry[k1]});
    }
  }
  mesh.refine(max_len);
  return mesh.to_list();
}

// What a renderer needs to run the forward projection and the quad layout
// itself. 'constants': tan, cos of the edge angle, cot 30 degrees, sin, cos
// and value of the vertex angle G, R', R'^2, the face-plane origin (x, y),
// the face edge, and the face projection (0 ISEA, 1 Fuller, 2 IVEA). 'faces': 16
// numbers per face: the centre as a unit vector, then the two unit vectors
// along which the face's azimuth is read (azimuth = atan2(p . b, p . a)), each
// followed by a zero, then the face's quad, its 60-degree turns into the quad
// and the offset after them. 'edges': 36 integers per quad: the quad of its
// far corner, whether it is a vertex quad, two zeros, then for each edge in
// QuadEdge order its QuadEdgeMap: the quad across it, the vertex quad of a
// far edge's starting corner (-1 for none), and the six coefficients.
// 'n_faces' and 'n_quads' count the solid's faces and quads.
// [[Rcpp::export]]
List cpp_globe_projection(NumericVector icosa) {
  activate_grid(icosa);
  const PolyData& S = poly();
  const SolidTopology& t = *S.topo;
  const SnyderParams& k = t.snyder;
  NumericVector constants = NumericVector::create(
      k.tan_el, k.cos_el, k.cot_30, k.sin_g, k.cos_g, k.g_angle,
      k.r1, k.r1_squared, k.origin_x_off, k.origin_y_off, k.edge,
      static_cast<double>(active_projection()));

  NumericVector faces(t.n_faces * 16);
  for (int f = 0; f < t.n_faces; f++) {
    const double lon = S.centers[f].lon, lat = S.centers[f].lat;
    const double centre[3] = {std::cos(lat) * std::cos(lon),
                              std::cos(lat) * std::sin(lon), std::sin(lat)};
    const double north[3] = {-std::sin(lat) * std::cos(lon),
                             -std::sin(lat) * std::sin(lon), std::cos(lat)};
    const double east[3] = {-std::sin(lon), std::cos(lon), 0.0};
    // project_core() reads the azimuth from north towards east and turns it
    // back by the face's offset.
    const double ca = std::cos(S.face_azimuth_offset[f]);
    const double sa = std::sin(S.face_azimuth_offset[f]);
    double* row = &faces[16 * f];
    for (int r = 0; r < 3; r++) {
      row[r] = centre[r];
      row[4 + r] = north[r] * ca + east[r] * sa;
      row[8 + r] = east[r] * ca - north[r] * sa;
    }
    const FacePlacement& pl = t.placement[f];
    row[12] = pl.quad;
    row[13] = pl.rotations;
    row[14] = pl.offset_x;
    row[15] = pl.offset_y;
  }

  IntegerVector edges(t.n_quads() * 36);
  for (int q = 0; q < t.n_quads(); q++) {
    int* row = &edges[36 * q];
    row[0] = t.corner[q][kCornerFar];
    row[1] = t.is_pole(q) ? 1 : 0;
    if (t.is_pole(q)) continue;
    for (int e = 0; e < 4; e++) {
      const QuadEdgeMap& m = t.edge[q][e];
      int* em = row + 4 + 8 * e;
      em[0] = m.quad;
      em[1] = m.pole;
      for (int c = 0; c < 2; c++) {
        for (int d = 0; d < 3; d++) em[2 + 3 * c + d] = m.k[c][d];
      }
    }
  }

  return List::create(_["constants"] = constants, _["faces"] = faces,
                      _["edges"] = edges, _["n_faces"] = t.n_faces,
                      _["n_quads"] = t.n_quads());
}

// Polygons of the sphere as one mesh on the faces. 'polygons' is a list of
// polygons, each a list of closed rings (n x 2 matrices of lon/lat, outer
// ring first, then holes). Each polygon is triangulated in the gnomonic
// projection about its centre, where its great-circle edges are straight;
// every triangle is cut into its parts on the sectors of the faces it crosses
// (sector_normals), and the parts are refined along great circles until no
// edge is longer than 'max_len'.
// Item k is polygon k.
// [[Rcpp::export]]
List cpp_globe_polygons(NumericVector icosa, List polygons, double max_len) {
  activate_icosa(icosa);
  FaceMesh mesh(/*great_circle=*/true);
  EarSlicer slicer;
  std::vector<double> gx, gy;
  std::vector<V3> pts, part, whole;
  std::vector<int> ring_start, tri;
  std::vector<double> tx, ty;

  for (R_xlen_t k = 0; k < polygons.size(); k++) {
    List rings = polygons[k];
    pts.clear();
    ring_start.clear();
    for (R_xlen_t r = 0; r < rings.size(); r++) {
      NumericMatrix ring = rings[r];
      int n = ring.nrow();
      // The ring is closed; its last point repeats the first.
      if (n >= 2 && ring(0, 0) == ring(n - 1, 0) && ring(0, 1) == ring(n - 1, 1)) n--;
      if (n < 3) {
        if (r == 0) break;
        continue;
      }
      ring_start.push_back(static_cast<int>(pts.size()));
      for (int i = 0; i < n; i++) pts.push_back(lonlat_vec(ring(i, 0), ring(i, 1)));
    }
    if (ring_start.empty()) continue;

    int n_outer = ring_start.size() > 1 ? ring_start[1] : static_cast<int>(pts.size());
    V3 c{0.0, 0.0, 0.0};
    for (int i = 0; i < n_outer; i++) c = c + pts[i];
    if (dot(c, c) < 1e-20) stop("polygon %d has no centre on the sphere", k + 1);
    c = unit(c);
    V3 u = unit(cross(std::fabs(c.z) < 0.9 ? V3{0, 0, 1} : V3{1, 0, 0}, c));
    V3 v = cross(c, u);
    gx.resize(pts.size());
    gy.resize(pts.size());
    for (size_t i = 0; i < pts.size(); i++) {
      double d = dot(pts[i], c);
      if (d < 0.02) {
        stop("polygon %d reaches more than 88 degrees from its centre; "
             "split it into smaller polygons", k + 1);
      }
      gx[i] = dot(pts[i], u) / d;
      gy[i] = dot(pts[i], v) / d;
    }

    tri.clear();
    slicer.run(gx, gy, ring_start, tri);

    for (size_t t = 0; t < tri.size(); t += 3) {
      V3 a = pts[tri[t]], b = pts[tri[t + 1]], d = pts[tri[t + 2]];
      for (int f = 0; f < poly().n_faces(); f++) {
        part.assign({a, b, d});
        clip_to_planes(face_normals(f), part);
        if (part.size() < 3) continue;
        whole.swap(part);
        for (int sector = 0; sector < n_face_sectors(); sector++) {
          part = whole;
          clip_to_planes(sector_normals(f, sector), part);
          if (part.size() < 3) continue;
          tx.clear();
          ty.clear();
          for (const V3& p : part) {
            auto t2 = face_tri(f, p);
            tx.push_back(t2.first);
            ty.push_back(t2.second);
          }
          int first = -1, prev = -1;
          for (size_t i = 0; i < part.size(); i++) {
            double s[3] = {part[i].x, part[i].y, part[i].z};
            int id = mesh.vertex(static_cast<int>(k + 1), f, tx[i], ty[i], s);
            if (i == 0) {
              first = id;
            } else if (i >= 2) {
              mesh.triangle(first, prev, id);
            }
            prev = id;
          }
        }
      }
    }
  }
  mesh.refine(max_len);
  return mesh.to_list();
}

// Polylines of the sphere placed on the faces. Points are lon/lat, 'path'
// numbers the polyline each point belongs to (consecutive points of one
// polyline). Segments are great-circle arcs, split where they cross a face
// edge and every 'max_angle' radians. Returns one row per point, with columns
// path, face (from 0), solid x/y/z, sphere x/y/z and triangle coordinates
// tx, ty, in the layout of cpp_cell_surface_paths(). A point where a path
// crosses from one face to the next appears on both faces, so the points of
// one face run unbroken up to the face edge, also in a layout of the unfolded
// solid that puts the two faces apart.
// [[Rcpp::export]]
NumericMatrix cpp_sphere_paths_on_faces(NumericVector icosa,
                                        NumericVector lon, NumericVector lat,
                                        IntegerVector path, double max_angle) {
  activate_icosa(icosa);
  if (!(max_angle > 0.0)) stop("max_angle must be positive");
  R_xlen_t n = lon.size();
  std::vector<double> rows;
  double solid[3];
  auto emit = [&](int id, int face, V3 p) {
    auto t = face_tri(face, p);
    face_tri_to_solid(face, t.first, t.second, solid);
    rows.insert(rows.end(), {static_cast<double>(id), static_cast<double>(face),
                             solid[0], solid[1], solid[2], p.x, p.y, p.z,
                             t.first, t.second});
  };

  for (R_xlen_t k = 0; k < n; k++) {
    V3 a = lonlat_vec(lon[k], lat[k]);
    bool last = k + 1 == n || path[k + 1] != path[k];
    if (last) {
      emit(path[k], sphere_face(a), a);
      continue;
    }
    V3 b = lonlat_vec(lon[k + 1], lat[k + 1]);
    double s = 0.0;
    int face = sphere_face(a);
    // A segment shorter than a half circle crosses few faces.
    for (int guard = 0; guard < 20; guard++) {
      const auto nrm = face_normals(face);
      double t_out = 1.0;
      for (int e = 0; e < 3; e++) {
        double da = dot(nrm[e], a), db = dot(nrm[e], b);
        if (db < 0.0 && da > db) {
          double t = da / (da - db);
          if (t > s - 1e-12 && t < t_out) t_out = t;
        }
      }
      if (t_out < s) t_out = s;
      V3 p0 = unit(a + s * (b - a)), p1 = unit(a + t_out * (b - a));
      double angle = std::acos(std::max(-1.0, std::min(1.0, dot(p0, p1))));
      int steps = std::max(1, static_cast<int>(std::ceil(angle / max_angle)));
      for (int i = 0; i < steps; i++) {
        double t = s + (t_out - s) * i / steps;
        emit(path[k], face, unit(a + t * (b - a)));
      }
      if (t_out >= 1.0) {
        // The next segment starts on another face when b lies on this
        // face's edge; b then also closes this face's run of points.
        if (sphere_face(b) != face) emit(path[k], face, b);
        break;
      }
      emit(path[k], face, p1);
      s = t_out;
      V3 next = unit(a + std::min(1.0, s + 1e-9) * (b - a));
      int best = -1;
      double best_m = -std::numeric_limits<double>::infinity();
      for (int f = 0; f < poly().n_faces(); f++) {
        if (f == face) continue;
        double m = face_margin(f, next);
        if (m > best_m) {
          best_m = m;
          best = f;
        }
      }
      face = best;
      if (s >= 1.0 - 1e-12) break;
    }
  }

  constexpr int n_col = 10;
  R_xlen_t n_row = static_cast<R_xlen_t>(rows.size() / n_col);
  NumericMatrix out(n_row, n_col);
  for (R_xlen_t r = 0; r < n_row; r++) {
    for (int col = 0; col < n_col; col++) out(r, col) = rows[r * n_col + col];
  }
  colnames(out) = CharacterVector::create("cell", "face", "solid_x", "solid_y",
                                          "solid_z", "sphere_x", "sphere_y",
                                          "sphere_z", "tx", "ty");
  return out;
}

static const char* BASE64_TABLE =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string hexify::base64_encode(const unsigned char* bytes, size_t n) {
  const char* table = BASE64_TABLE;
  std::string out;
  out.reserve((n + 2) / 3 * 4);
  size_t i = 0;
  for (; i + 2 < n; i += 3) {
    uint32_t v = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
    out += table[(v >> 18) & 63];
    out += table[(v >> 12) & 63];
    out += table[(v >> 6) & 63];
    out += table[v & 63];
  }
  if (i < n) {
    uint32_t v = bytes[i] << 16;
    if (i + 1 < n) v |= bytes[i + 1] << 8;
    out += table[(v >> 18) & 63];
    out += table[(v >> 12) & 63];
    out += i + 1 < n ? table[(v >> 6) & 63] : '=';
    out += '=';
  }
  return out;
}

std::string hexify::base64_words(const std::vector<uint32_t>& words) {
  std::vector<unsigned char> bytes(4 * words.size());
  for (size_t k = 0; k < words.size(); k++) {
    for (int b = 0; b < 4; b++) bytes[4 * k + b] = (words[k] >> (8 * b)) & 0xFF;
  }
  return base64_encode(bytes.data(), bytes.size());
}

// Numbers as the base64 text of their little-endian bytes: 32-bit floats for
// type "f32", unsigned 32-bit integers for "u32".
// [[Rcpp::export]]
std::string cpp_base64_buffer(NumericVector x, std::string type) {
  bool f32 = type == "f32";
  if (!f32 && type != "u32") stop("type must be \"f32\" or \"u32\"");
  R_xlen_t n = x.size();
  std::vector<uint32_t> words(static_cast<size_t>(n));
  for (R_xlen_t k = 0; k < n; k++) {
    if (f32) {
      float v = static_cast<float>(x[k]);
      std::memcpy(&words[k], &v, 4);
    } else {
      words[k] = static_cast<uint32_t>(x[k]);
    }
  }
  return base64_words(words);
}

// Bytes as base64 text
// [[Rcpp::export]]
std::string cpp_base64_bytes(RawVector x) {
  return base64_encode(x.begin(), static_cast<size_t>(x.size()));
}

// The bytes of base64 text; padding and characters outside the alphabet end
// or are skipped.
// [[Rcpp::export]]
RawVector cpp_base64_decode(std::string text) {
  int value[256];
  std::fill(value, value + 256, -1);
  for (int k = 0; k < 64; k++) value[static_cast<unsigned char>(BASE64_TABLE[k])] = k;
  std::vector<unsigned char> out;
  out.reserve(text.size() / 4 * 3);
  uint32_t acc = 0;
  int bits = 0;
  for (unsigned char c : text) {
    if (c == '=') break;
    int v = value[c];
    if (v < 0) continue;
    acc = (acc << 6) | static_cast<uint32_t>(v);
    bits += 6;
    if (bits >= 8) {
      bits -= 8;
      out.push_back(static_cast<unsigned char>((acc >> bits) & 0xFF));
    }
  }
  RawVector r(out.size());
  std::copy(out.begin(), out.end(), r.begin());
  return r;
}
