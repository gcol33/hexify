#ifndef HEX_SEQNUM_H
#define HEX_SEQNUM_H

#include <cstdint>

namespace hexify {

// Sequential numbering: converts between (face, i, j, resolution) and integer seqnum
// This provides PLANE addressing mode compatible with DGGRID

// ============================================================================
// APERTURE 3
// ============================================================================

// Convert (face, i, j, resolution) → integer seqnum
uint64_t cell_to_seqnum_ap3(int face, long long i, long long j, int resolution);

// Convert integer seqnum → (face, i, j) at given resolution
void seqnum_to_cell_ap3(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j);

// ============================================================================
// APERTURE 4
// ============================================================================

uint64_t cell_to_seqnum_ap4(int face, long long i, long long j, int resolution);

void seqnum_to_cell_ap4(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j);

// ============================================================================
// APERTURE 7
// ============================================================================

uint64_t cell_to_seqnum_ap7(int face, long long i, long long j, int resolution);

void seqnum_to_cell_ap7(uint64_t seqnum, int resolution,
                        int& face, long long& i, long long& j);

} // namespace hexify

#endif // HEX_SEQNUM_H
