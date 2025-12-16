// cell_numbering.h
// Cell ID Assignment for ISEA Hexagonal Grids
//
// Converts between (quad, i, j) grid coordinates and sequential cell IDs.
// Cell IDs provide a compact 1D representation of cells for storage/indexing.
//
// The cell ID scheme is a simple bijective mapping:
//   cell_id = quad * cells_per_quad + local_index + 1
//
// This is a standard approach for discrete global grids documented in:
//   Sahr, K., White, D., & Kimerling, A.J. (2003). "Geodesic Discrete
//   Global Grid Systems." Cartography and Geographic Information Science.
//
// Cell count formula: N = 10 * aperture^resolution + 2
// The "+2" accounts for the two pentagon cells at poles.
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#ifndef HEXIFY_CELL_NUMBERING_H
#define HEXIFY_CELL_NUMBERING_H

#include <cstdint>

namespace hexify {

// ============================================================================
// Unified API
// ============================================================================

/**
 * Convert (quad, i, j) coordinates to cell ID.
 *
 * @param quad       Quad index (0-11)
 * @param i          Column index within quad
 * @param j          Row index within quad
 * @param aperture   Grid aperture (3, 4, or 7)
 * @param resolution Grid resolution level
 * @return           1-based cell ID
 */
uint64_t quad_ij_to_cell_id(int quad, long long i, long long j,
                            int aperture, int resolution);

/**
 * Convert cell ID to (quad, i, j) coordinates.
 *
 * @param cell_id    1-based cell ID
 * @param aperture   Grid aperture (3, 4, or 7)
 * @param resolution Grid resolution level
 * @param[out] quad  Quad index (0-11)
 * @param[out] i     Column index within quad
 * @param[out] j     Row index within quad
 */
void cell_id_to_quad_ij(uint64_t cell_id, int aperture, int resolution,
                        int& quad, long long& i, long long& j);

// ============================================================================
// Aperture-Specific Functions (for performance-critical code paths)
// ============================================================================

// Aperture 3
uint64_t quad_ij_to_cell_id_ap3(int quad, long long i, long long j, int resolution);
void cell_id_to_quad_ij_ap3(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j);

// Aperture 4
uint64_t quad_ij_to_cell_id_ap4(int quad, long long i, long long j, int resolution);
void cell_id_to_quad_ij_ap4(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j);

// Aperture 7
uint64_t quad_ij_to_cell_id_ap7(int quad, long long i, long long j, int resolution);
void cell_id_to_quad_ij_ap7(uint64_t cell_id, int resolution,
                            int& quad, long long& i, long long& j);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Get total number of cells at given resolution and aperture.
 * Formula: N = 10 * aperture^resolution + 2
 */
uint64_t total_cell_count(int aperture, int resolution);

/**
 * Get maximum valid cell ID at given resolution and aperture.
 */
uint64_t max_cell_id(int aperture, int resolution);

/**
 * Validate that a cell ID is within valid range.
 */
bool is_valid_cell_id(uint64_t cell_id, int aperture, int resolution);

} // namespace hexify

#endif // HEXIFY_CELL_NUMBERING_H
