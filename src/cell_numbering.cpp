// cell_numbering.cpp
// Cell ID Assignment for ISEA Hexagonal Grids
//
// Implements bijective mapping between (quad, i, j) coordinates and cell IDs.
//
// The mapping uses a simple linearization scheme:
//   cell_id = quad * cells_per_quad + row_major_index(i, j) + 1
//
// Cell IDs are 1-based to match common GIS conventions and allow 0 to
// represent "no cell" or "invalid" in application code.
//
// References:
//   Sahr, K., White, D., & Kimerling, A.J. (2003). "Geodesic Discrete
//   Global Grid Systems." Cartography and Geographic Information Science.
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "cell_numbering.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// ============================================================================
// Internal Helper Functions
// ============================================================================

/**
 * Compute the grid dimension (max coordinate + 1) for a given aperture/resolution.
 *
 * The grid dimension determines the range of valid (i, j) coordinates.
 * For each aperture, the dimension grows as a power of the scale factor:
 *   - Aperture 3: dimension = 3^((res+1)/2) due to Class I/II alternation
 *   - Aperture 4: dimension = 2^res
 *   - Aperture 7: dimension = 7^res
 */
inline long long grid_dimension(int aperture, int resolution) {
    if (resolution <= 0) return 1;

    long long dim = 1;
    switch (aperture) {
        case 3: {
            // Aperture 3 uses effective resolution (res+1)/2 due to
            // alternating Class I/II quantization
            int eff_res = (resolution + 1) / 2;
            for (int r = 0; r < eff_res; ++r) dim *= 3;
            break;
        }
        case 4:
            for (int r = 0; r < resolution; ++r) dim *= 2;
            break;
        case 7:
            for (int r = 0; r < resolution; ++r) dim *= 7;
            break;
        default:
            break;
    }
    return dim;
}

/**
 * Compute cells per quad at given aperture/resolution.
 * Each quad contains dim * dim cells in a row-major layout.
 */
inline uint64_t cells_per_quad(int aperture, int resolution) {
    long long dim = grid_dimension(aperture, resolution);
    return static_cast<uint64_t>(dim) * static_cast<uint64_t>(dim);
}

/**
 * Convert (i, j) to linear index within a quad using row-major order.
 * index = i * dimension + j
 */
inline uint64_t ij_to_linear(long long i, long long j, long long dim) {
    // Clamp to valid range
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= dim) i = dim - 1;
    if (j >= dim) j = dim - 1;

    return static_cast<uint64_t>(i) * static_cast<uint64_t>(dim) +
           static_cast<uint64_t>(j);
}

/**
 * Convert linear index back to (i, j).
 */
inline void linear_to_ij(uint64_t index, long long dim,
                         long long& i, long long& j) {
    i = static_cast<long long>(index / static_cast<uint64_t>(dim));
    j = static_cast<long long>(index % static_cast<uint64_t>(dim));
}

/**
 * Core implementation of quad_ij -> cell_id conversion.
 */
inline uint64_t quad_ij_to_cell_id_impl(int quad, long long i, long long j,
                                        int aperture, int resolution) {
    // Resolution 0: just the quad number + 1 (12 quads)
    if (resolution == 0) {
        return static_cast<uint64_t>(quad) + 1;
    }

    long long dim = grid_dimension(aperture, resolution);
    uint64_t cpq = static_cast<uint64_t>(dim) * static_cast<uint64_t>(dim);
    uint64_t local_idx = ij_to_linear(i, j, dim);

    // cell_id = quad * cells_per_quad + local_index + 1
    return static_cast<uint64_t>(quad) * cpq + local_idx + 1;
}

/**
 * Core implementation of cell_id -> quad_ij conversion.
 */
inline void cell_id_to_quad_ij_impl(uint64_t cell_id, int aperture, int resolution,
                                    int& quad, long long& i, long long& j) {
    // Resolution 0: quad = cell_id - 1
    if (resolution == 0) {
        quad = static_cast<int>(cell_id - 1);
        i = 0;
        j = 0;
        return;
    }

    // Convert to 0-based index
    uint64_t idx = cell_id - 1;

    long long dim = grid_dimension(aperture, resolution);
    uint64_t cpq = static_cast<uint64_t>(dim) * static_cast<uint64_t>(dim);

    quad = static_cast<int>(idx / cpq);
    uint64_t local_idx = idx % cpq;

    linear_to_ij(local_idx, dim, i, j);
}

} // anonymous namespace

// ============================================================================
// Unified API Implementation
// ============================================================================

uint64_t quad_ij_to_cell_id(int quad, long long i, long long j,
                            int aperture, int resolution) {
    if (quad < 0 || quad > 11) {
        throw std::runtime_error("quad_ij_to_cell_id: quad must be 0-11");
    }
    if (resolution < 0) {
        throw std::runtime_error("quad_ij_to_cell_id: resolution must be >= 0");
    }
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        throw std::runtime_error("quad_ij_to_cell_id: aperture must be 3, 4, or 7");
    }

    return quad_ij_to_cell_id_impl(quad, i, j, aperture, resolution);
}

void cell_id_to_quad_ij(uint64_t cell_id, int aperture, int resolution,
                        int& quad, long long& i, long long& j) {
    if (cell_id < 1) {
        throw std::runtime_error("cell_id_to_quad_ij: cell_id must be >= 1");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_id_to_quad_ij: resolution must be >= 0");
    }
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        throw std::runtime_error("cell_id_to_quad_ij: aperture must be 3, 4, or 7");
    }

    cell_id_to_quad_ij_impl(cell_id, aperture, resolution, quad, i, j);

    if (quad > 11) {
        throw std::runtime_error("cell_id_to_quad_ij: cell_id out of range");
    }
}

// ============================================================================
// Aperture-Specific Functions
// ============================================================================
//
// These provide slightly faster paths by avoiding the aperture switch,
// useful for hot loops processing millions of cells.

// --- Aperture 3 ---

uint64_t quad_ij_to_cell_id_ap3(int quad, long long i, long long j, int resolution) {
    if (quad < 0 || quad > 11) {
        throw std::runtime_error("quad_ij_to_cell_id_ap3: quad must be 0-11");
    }
    if (resolution < 0) {
        throw std::runtime_error("quad_ij_to_cell_id_ap3: resolution must be >= 0");
    }
    return quad_ij_to_cell_id_impl(quad, i, j, 3, resolution);
}

void cell_id_to_quad_ij_ap3(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j) {
    if (cell_id < 1) {
        throw std::runtime_error("cell_id_to_quad_ij_ap3: cell_id must be >= 1");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_id_to_quad_ij_ap3: resolution must be >= 0");
    }

    cell_id_to_quad_ij_impl(cell_id, 3, resolution, quad, i, j);

    if (quad > 11) {
        throw std::runtime_error("cell_id_to_quad_ij_ap3: cell_id out of range");
    }
}

// --- Aperture 4 ---

uint64_t quad_ij_to_cell_id_ap4(int quad, long long i, long long j, int resolution) {
    if (quad < 0 || quad > 11) {
        throw std::runtime_error("quad_ij_to_cell_id_ap4: quad must be 0-11");
    }
    if (resolution < 0) {
        throw std::runtime_error("quad_ij_to_cell_id_ap4: resolution must be >= 0");
    }
    return quad_ij_to_cell_id_impl(quad, i, j, 4, resolution);
}

void cell_id_to_quad_ij_ap4(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j) {
    if (cell_id < 1) {
        throw std::runtime_error("cell_id_to_quad_ij_ap4: cell_id must be >= 1");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_id_to_quad_ij_ap4: resolution must be >= 0");
    }

    cell_id_to_quad_ij_impl(cell_id, 4, resolution, quad, i, j);

    if (quad > 11) {
        throw std::runtime_error("cell_id_to_quad_ij_ap4: cell_id out of range");
    }
}

// --- Aperture 7 ---

uint64_t quad_ij_to_cell_id_ap7(int quad, long long i, long long j, int resolution) {
    if (quad < 0 || quad > 11) {
        throw std::runtime_error("quad_ij_to_cell_id_ap7: quad must be 0-11");
    }
    if (resolution < 0) {
        throw std::runtime_error("quad_ij_to_cell_id_ap7: resolution must be >= 0");
    }
    return quad_ij_to_cell_id_impl(quad, i, j, 7, resolution);
}

void cell_id_to_quad_ij_ap7(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j) {
    if (cell_id < 1) {
        throw std::runtime_error("cell_id_to_quad_ij_ap7: cell_id must be >= 1");
    }
    if (resolution < 0) {
        throw std::runtime_error("cell_id_to_quad_ij_ap7: resolution must be >= 0");
    }

    cell_id_to_quad_ij_impl(cell_id, 7, resolution, quad, i, j);

    if (quad > 11) {
        throw std::runtime_error("cell_id_to_quad_ij_ap7: cell_id out of range");
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

uint64_t total_cell_count(int aperture, int resolution) {
    if (resolution < 0) return 0;
    if (resolution == 0) return 12;  // 12 quads at resolution 0

    // N = 12 * cells_per_quad
    return 12 * cells_per_quad(aperture, resolution);
}

uint64_t max_cell_id(int aperture, int resolution) {
    return total_cell_count(aperture, resolution);
}

bool is_valid_cell_id(uint64_t cell_id, int aperture, int resolution) {
    if (cell_id < 1) return false;
    return cell_id <= max_cell_id(aperture, resolution);
}

} // namespace hexify
