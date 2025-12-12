// cell_numbering.h
// Cell ID Assignment for ISEA Hexagonal Grids
//
// Converts between (face, i, j) grid coordinates and sequential cell IDs.
// Cell IDs provide a compact 1D representation of cells for storage/indexing.
//
// Formula: Total cells = 10 * aperture^resolution + 2
// The "+2" accounts for the two pentagon cells at poles.
//
// References:
// - Sahr et al. (2003) "Geodesic Discrete Global Grid Systems"
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#ifndef HEXIFY_CELL_NUMBERING_H
#define HEXIFY_CELL_NUMBERING_H

#include <cstdint>

namespace hexify {

// ============================================================================
// APERTURE 3
// ============================================================================

// Convert (face, i, j, resolution) → cell ID
uint64_t cell_to_seqnum_ap3(int face, long long i, long long j, int resolution);

// Convert cell ID → (face, i, j) at given resolution
void seqnum_to_cell_ap3(uint64_t cell_id, int resolution,
                        int& face, long long& i, long long& j);

// ============================================================================
// APERTURE 4
// ============================================================================

uint64_t cell_to_seqnum_ap4(int face, long long i, long long j, int resolution);

void seqnum_to_cell_ap4(uint64_t cell_id, int resolution,
                        int& face, long long& i, long long& j);

// ============================================================================
// APERTURE 7
// ============================================================================

uint64_t cell_to_seqnum_ap7(int face, long long i, long long j, int resolution);

void seqnum_to_cell_ap7(uint64_t cell_id, int resolution,
                        int& face, long long& i, long long& j);

} // namespace hexify

#endif // HEXIFY_CELL_NUMBERING_H
