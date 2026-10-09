#pragma once
// Hex9 cells on hexify's octahedron: between a cell's lattice point (hex9.h)
// and the face triangle coordinates and quad lattice coordinates the rest of
// the package reads.
//
// The octahedron's standard frame puts vertex 0 on +Z, vertices 1 and 2 on +X
// and +Y, vertices 3 and 4 on -X and -Y and vertex 5 on -Z (polyhedron.cpp),
// so the Hex9 lattice is fixed to the solid and turns with its orientation.
// On a face, a lattice point's weights on the face's three vertices are its
// absolute coordinates on their axes; in the face's triangle coordinates the
// vertices sit at the top, lower left and lower right.
//
// A level-L cell is stored in the quads as the aperture-3 Class II lattice of
// resolution 2L + 1 stores its cells, by the substrate coordinates of its
// centre on a quad edge of 3^(L+1) steps; the Hex9 centres are a coset of
// that lattice, off the vertices of the solid.

#include "hex9.h"

namespace hexify {
namespace hex9 {

// The aperture-3 resolution whose quad frame stores level `level`
inline int frame_resolution(int level) { return 2 * level + 1; }

// The cell of `level` holding the point (tx, ty) of `face` of the active
// solid, which must be the octahedron.
OctPoint face_point_cell(int face, double tx, double ty, int level);

// A lattice point of `level` as a face holding it and its triangle
// coordinates there.
void lattice_face_xy(const OctPoint& p, int level, int& face, double& tx, double& ty);

// A cell's centre as the quad that owns it and its substrate (i, j) there
void cell_quad_ij(const OctPoint& centre, int level, int& quad, long long& i, long long& j);

// The lattice point at substrate (i, j) of `quad`, a point inside the quad's
// box. False when (i, j) is no lattice point of a face.
bool quad_ij_lattice(int quad, long long i, long long j, int level, OctPoint& p);

} // namespace hex9
} // namespace hexify
