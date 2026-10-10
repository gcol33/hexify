#pragma once
// The trained area-correcting warp of Hex9's own projection (Griffin 2026,
// Sec. 11b): libhex9's Sphere-L6 field, a displacement of each octant face
// stored at the vertices of the face's level-6 lattice (3^7 steps along an
// edge) with its gradients, on the sixth of the face its symmetries leave
// (the wedge), and read between vertices by a Clough-Tocher interpolant.
//
// With Kaseorg's projection AK, a point on the sphere goes to the face point
// P = AK^-1(p), and the grid's lattice point L under it is the solution of
// L + d(L) = P; back, a lattice point L goes to the sphere at AK(L + d(L)).
//
// The field comes as libhex9's .h9warp file (version 4, the fundamental
// wedge with gradients), read here byte for byte; it is not shipped with
// hexify.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hexify {
namespace hex9 {

// Load the field from the bytes of a .h9warp file. False, with `err`, when
// the bytes are not a version-4 wedge field with gradients or fail its CRC.
bool warp_load(const unsigned char* data, std::size_t len, std::string& err);

// Whether a field is loaded
bool warp_ready();

// The displacement d at face-plane point (x, y) of any octant face, in the
// face's triangle coordinates (unit edge, origin at the lower-left vertex).
// The field has the face's six symmetries, so it reads the same on every
// face whichever vertex the face's coordinates start from.
void warp_delta(double x, double y, double& dx, double& dy);

// L + d(L): the lattice point (x, y) to the point AK maps to the sphere
void warp_apply(double x, double y, double& px, double& py);

// The lattice point L with L + d(L) = (px, py), by Newton's method
void warp_solve(double px, double py, double& x, double& y);

// The field as the globe shader reads it (globe.wgsl's warp_word()): n, the
// first point of each wedge column, then each point's six numbers as 32-bit
// floats; empty when no field is loaded
std::vector<uint32_t> warp_texture_words();

} // namespace hex9
} // namespace hexify
