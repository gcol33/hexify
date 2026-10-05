// globe_mesh.h
// Triangle meshes and polylines on the icosahedron and the sphere at once
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once
#include <Rcpp.h>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace hexify {

// A triangle mesh laid on the faces of the icosahedron. Every vertex sits on
// one face, in that face's triangle coordinates, so it has a place on the
// flat face and one on the sphere. Each triangle lies within one face and
// carries an item number (a cell or a polygon), shared by its vertices.
class FaceMesh {
public:
  // With 'great_circle', an edge is split at the midpoint of its great-circle
  // arc on the sphere, and every vertex must be given its sphere position;
  // otherwise at the midpoint of the straight line on its face.
  explicit FaceMesh(bool great_circle = false) : great_circle_(great_circle) {}

  // The vertex (tx, ty) of a face for an item; the same point of the same
  // face and item is stored once. 'sphere', when given, is the vertex's place
  // on the unit sphere, kept instead of the inverse projection.
  int vertex(int item, int face, double tx, double ty,
             const double* sphere = nullptr);
  void triangle(int a, int b, int c);

  // A convex polygon of one face, in its triangle coordinates, as a fan.
  void convex_polygon(int item, int face, const std::vector<double>& tx,
                      const std::vector<double>& ty);

  // Split every edge longer than 'max_len' (in triangle coordinates, a face
  // edge is 1) at its midpoint until none is. The split of an edge depends on
  // the edge alone, so two triangles that share an edge split it alike and
  // the mesh keeps no gaps.
  void refine(double max_len);

  // The mesh for R: solid and sphere positions (x, y, z per vertex), the
  // vertex indices of each triangle (from 0) and the item of each vertex.
  Rcpp::List to_list() const;

  size_t n_vertices() const { return item_.size(); }

private:
  struct Key {
    int item, face;
    double tx, ty;
    bool operator==(const Key& o) const {
      return item == o.item && face == o.face && tx == o.tx && ty == o.ty;
    }
  };
  struct KeyHash {
    size_t operator()(const Key& k) const;
  };

  bool great_circle_;
  std::vector<int> item_, face_;
  std::vector<double> tx_, ty_;
  std::vector<double> sphere_;      // x, y, z per vertex; NaN when not given
  std::vector<int> tri_;
  std::unordered_map<Key, int, KeyHash> index_;
};

// The corners of a face in its triangle coordinates, in the order of the
// face's vertices.
void face_tri_corners(int face, double tx[3], double ty[3]);

// The point (tx, ty) of face 'from' in the triangle coordinates of face
// 'to', read through the projection of 'to' extended past its edges.
void face_tri_to_face(int from, double tx, double ty, int to,
                      double& out_x, double& out_y);

// Cut a polygon given in a face's triangle coordinates to that triangle
// (Sutherland-Hodgman clipping against its three edges).
void clip_to_face_tri(int face, std::vector<double>& x, std::vector<double>& y);

} // namespace hexify
