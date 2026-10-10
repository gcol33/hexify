// globe_table.h
// The cells of a grid as the globe shader finds and reads them
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once
#include <cstdint>
#include <vector>
#include "polyhedron.h"

namespace hexify {

// One resolution of a grid as the globe shader reads it. Scaled by `dim`,
// the quad's side in substrate steps, a quad-plane point's nearest cell
// centre is its nearest multiple of the generator ga + gb*omega (omega =
// exp(2*pi*i/3)) in the substrate's (i, j); those multiples are the points
// with j = c * i (mod index). A quad holds `per_quad` cells, numbered as
// cell_index_2d() numbers that sublattice, and the grid `n_cells`, the two
// vertex quads' cells first and last. Aperture 7 stores surrogates but
// numbers its cells by their substrate centres, which is the same count.
//
// Hex9 (aperture 9, hex9.h) has no cell at a vertex of the solid: its cells
// are the points of a coset of that sublattice, j = c * i + coset[q] (mod
// index) in quad q, numbered quad by quad from the first diamond quad's.
// Its cell IDs are addresses, which the table reads through hex9_level.
struct GlobeFrame {
    long long dim;
    long long index;
    long long c;
    long long ga;
    long long gb;
    uint64_t per_quad;
    uint64_t n_cells;
    int hex9_level = -1;          // the Hex9 level, -1 for every other grid
    int coset[kMaxVerts] = {0};   // per quad; 0 but on Hex9
    bool vertex_cells() const { return hex9_level < 0; }
};

// The frame of a pure aperture at a resolution (empty `ap_seq`), or of a
// mixed sequence, whose resolution is one less than its length, on the
// active solid.
GlobeFrame globe_frame(int resolution, int aperture, const std::vector<int>& ap_seq);

} // namespace hexify
