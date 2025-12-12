// cell_numbering.cpp
// Cell ID Assignment for ISEA Hexagonal Grids
//
// Implements bijective mapping between (face, i, j) coordinates and cell IDs.
// Cell IDs are 1-based integers providing compact cell identification.
//
// Formula: face * cells_per_face + cell_index + 1
// where cell_index is computed from (i, j) as a linear index within the face.
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "cell_numbering.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

// ============================================================================
// Helper Functions
// ============================================================================

// Maximum coordinate value at given resolution for given aperture
// For aperture 3: max = 3^res - 1 (but grid is not square)
// For aperture 4: max = 2^res - 1
// For aperture 7: max varies based on resolution
static long long max_coord(int resolution, int aperture) {
    if (resolution == 0) return 0;
    long long m = 1;
    if (aperture == 3) {
        // For Z3, effRes = (res+1)/2
        int effRes = (resolution + 1) / 2;
        for (int r = 0; r < effRes; r++) m *= 3;
    } else if (aperture == 4) {
        for (int r = 0; r < resolution; r++) m *= 2;
    } else if (aperture == 7) {
        for (int r = 0; r < resolution; r++) m *= 7;
    }
    return m - 1;
}

// Cells per face at given resolution
// This is the number of valid cells in one face
static uint64_t cells_per_face(int resolution, int aperture) {
    if (resolution == 0) return 1;
    uint64_t n = 1;
    for (int r = 0; r < resolution; r++) {
        n *= aperture;
    }
    return n;
}

// Total cells across all faces
// For hex grids: 12 * cells_per_face (since we have 12 faces/quads)
static uint64_t total_cells(int resolution, int aperture) {
    if (resolution == 0) return 12;
    return 12 * cells_per_face(resolution, aperture);
}

// Convert (i, j) to a linear index within a face
// Uses row-major order: index = i * num_cols + j
// For Z3/Z7, the grid is not square, so we use the larger dimension
static uint64_t cell_to_linear_index(long long i, long long j, int resolution, int aperture) {
    if (resolution == 0) return 0;

    long long max_val = max_coord(resolution, aperture) + 1;

    // Clamp to valid range
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > max_val) i = max_val;
    if (j > max_val) j = max_val;

    return static_cast<uint64_t>(i) * max_val + j;
}

// Convert linear index back to (i, j)
static void linear_index_to_cell(uint64_t index, int resolution, int aperture,
                                  long long& i, long long& j) {
    if (resolution == 0) {
        i = 0;
        j = 0;
        return;
    }

    long long max_val = max_coord(resolution, aperture) + 1;
    i = index / max_val;
    j = index % max_val;
}

// ============================================================================
// Public API - Aperture 3
// ============================================================================

uint64_t cell_to_seqnum_ap3(int face, long long i, long long j, int resolution) {
    if (face < 0 || face > 19) {
        throw std::runtime_error("cell_numbering: face must be 0-19");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }

    // Resolution 0: just the face + 1 (20 faces)
    if (resolution == 0) {
        return static_cast<uint64_t>(face) + 1;
    }

    // Convert (i, j) to linear index within face
    uint64_t cell_idx = cell_to_linear_index(i, j, resolution, 3);

    // Total cells per face (using 20 triangular faces)
    long long max_val = max_coord(resolution, 3) + 1;
    uint64_t cpf = max_val * max_val;

    // Seqnum = face * cpf + cell_idx + 1 (1-based)
    return static_cast<uint64_t>(face) * cpf + cell_idx + 1;
}

void seqnum_to_cell_ap3(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j) {
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }
    if (seqnum < 1) {
        throw std::runtime_error("cell_numbering: seqnum must be >= 1");
    }

    // Resolution 0: face is seqnum - 1 (20 faces)
    if (resolution == 0) {
        face = static_cast<int>(seqnum - 1);
        if (face > 19) face = 19;
        i = 0;
        j = 0;
        return;
    }

    // Convert to 0-based
    uint64_t snum = seqnum - 1;

    // Cells per face (using 20 triangular faces)
    long long max_val = max_coord(resolution, 3) + 1;
    uint64_t cpf = max_val * max_val;

    // Calculate face and cell index
    face = static_cast<int>(snum / cpf);
    if (face > 19) face = 19;
    uint64_t cell_idx = snum % cpf;

    // Convert cell index to (i, j)
    linear_index_to_cell(cell_idx, resolution, 3, i, j);
}

// ============================================================================
// Public API - Aperture 4
// ============================================================================

uint64_t cell_to_seqnum_ap4(int face, long long i, long long j, int resolution) {
    if (face < 0 || face > 19) {
        throw std::runtime_error("cell_numbering: face must be 0-19");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }

    if (resolution == 0) {
        return static_cast<uint64_t>(face) + 1;
    }

    uint64_t cell_idx = cell_to_linear_index(i, j, resolution, 4);
    long long max_val = max_coord(resolution, 4) + 1;
    uint64_t cpf = max_val * max_val;

    return static_cast<uint64_t>(face) * cpf + cell_idx + 1;
}

void seqnum_to_cell_ap4(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j) {
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }
    if (seqnum < 1) {
        throw std::runtime_error("cell_numbering: seqnum must be >= 1");
    }

    if (resolution == 0) {
        face = static_cast<int>(seqnum - 1);
        if (face > 19) face = 19;
        i = 0;
        j = 0;
        return;
    }

    uint64_t snum = seqnum - 1;
    long long max_val = max_coord(resolution, 4) + 1;
    uint64_t cpf = max_val * max_val;

    face = static_cast<int>(snum / cpf);
    if (face > 19) face = 19;
    uint64_t cell_idx = snum % cpf;

    linear_index_to_cell(cell_idx, resolution, 4, i, j);
}

// ============================================================================
// Public API - Aperture 7
// ============================================================================

uint64_t cell_to_seqnum_ap7(int face, long long i, long long j, int resolution) {
    if (face < 0 || face > 19) {
        throw std::runtime_error("cell_numbering: face must be 0-19");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }

    if (resolution == 0) {
        return static_cast<uint64_t>(face) + 1;
    }

    uint64_t cell_idx = cell_to_linear_index(i, j, resolution, 7);
    long long max_val = max_coord(resolution, 7) + 1;
    uint64_t cpf = max_val * max_val;

    return static_cast<uint64_t>(face) * cpf + cell_idx + 1;
}

void seqnum_to_cell_ap7(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j) {
    if (resolution < 0) {
        throw std::runtime_error("cell_numbering: resolution must be >= 0");
    }
    if (seqnum < 1) {
        throw std::runtime_error("cell_numbering: seqnum must be >= 1");
    }

    if (resolution == 0) {
        face = static_cast<int>(seqnum - 1);
        if (face > 19) face = 19;
        i = 0;
        j = 0;
        return;
    }

    uint64_t snum = seqnum - 1;
    long long max_val = max_coord(resolution, 7) + 1;
    uint64_t cpf = max_val * max_val;

    face = static_cast<int>(snum / cpf);
    if (face > 19) face = 19;
    uint64_t cell_idx = snum % cpf;

    linear_index_to_cell(cell_idx, resolution, 7, i, j);
}

} // namespace hexify
