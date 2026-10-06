// index_z7.cpp
// Z7 Hierarchical Index Implementation for Aperture-7 Hexagonal Grids
//
// Encodes (quad, i, j) coordinates as a hierarchical digit string.
// The encoding traverses the aperture-7 hierarchy from coarse to fine,
// recording the child position (0-6) at each level.
//
// Algorithm:
// 1. Start at finest resolution, repeatedly coarsen using upAp7/upAp7r
// 2. At each level, compute child position as difference from parent center
// 3. Encode digits 0-6 representing the 7 child positions, led by the quad
//    and the unit digit the walk arrives at
//
// References:
// - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
// - H3 coordijk.c (Apache 2.0) for aperture-7 coordinate math
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "index_z7.h"
#include "polyhedron.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <vector>
#include <algorithm>

namespace hexify {
namespace z7 {

// ============================================================================
// Bijective aperture-7 hierarchical index (hexify-native)
// ============================================================================
// DGGRID's DgZ7StringRF encode is non-injective near the pentagon base cells:
// two geographically distinct cells can encode to one string (verified
// reproducer: lon/lat (5,45) and (-34.9,60.2) both -> "0045310"), so a
// bijective cell<->index round-trip is impossible with the faithful algorithm.
//
// These functions keep hexify's geographic quad fixed -- they reuse the exact
// aperture-7 digit machinery (upAp7/downAp7 + diffVec) but drop DGGRID's
// base-cell adjacency reassignment and pentagon digit-skip, the steps that
// merge distinct cells. Every (quad, i, j) then maps to a unique string and
// back. The string equals the DGGRID Z7 string for cells DGGRID does not
// reassign, and deviates only for the pentagon-region cells where DGGRID's own
// encoder collides. Input/output (i,j) are the Class I substrate coordinate.

// The level-0 coordinate a walk arrives at. A cell whose whole ancestry lies
// inside its quad arrives at the origin, and decoding from the origin recovers
// it. The quad is a rhombus while the aperture-7 parents are hexagons, so the
// quad boundary cuts through the parents of the cells along it; those arrive at
// one of the six neighbours of the origin instead, which the walk alone does
// not record. The arrival point is a unit digit, so the index carries it in its
// leading field as quad + n * digit, n the solid's number of quads (12 on the
// icosahedron): two digits still, and the plain quad DGGRID writes whenever
// the cell does not spill.
static IVec3D z7_seed_coord(int digit) {
    IVec3D ijk(0, 0, 0);
    ijk.neighbor((IVec3D::Direction) digit);
    return ijk;
}

// The base cell a seed names.
//
// The seed is the lattice point the hierarchy walk arrives at, at resolution 0.
// The quad's own origin means the whole ancestry lies inside the quad; any
// other arrival has crossed the quad's rhombic boundary towards one of its
// other three corners, whose vertex is the base cell DGGRID's Z7 encoder
// reassigns such a cell to. A vertex quad has no other corners.
static int seed_base_cell(int quadNum, int seed) {
    const SolidTopology& t = topo();
    if (t.is_pole(quadNum)) return quadNum;
    const IVec3D arrival = z7_seed_coord(seed);
    if (arrival.i() == 1) {
        return (arrival.j() == 0) ? t.corner[quadNum][kCornerI]
                                  : t.corner[quadNum][kCornerFar];
    }
    if (arrival.j() == 1) {
        return t.corner[quadNum][kCornerJ];
    }
    return quadNum;
}

// Walk the aperture-7 hierarchy from a Class I substrate coordinate up to
// resolution 0, recording the child ordinal of each level in digits[1..res].
// The return value is the coordinate the walk arrives at: the origin when every
// ancestor of the cell lies in the same quad, and a neighbouring lattice point
// when the quad's rhombic boundary cuts through one of them.
static IVec3D z7_walk_up(long long i, long long j, int resolution,
                         std::vector<IVec3D::Direction>& digits) {
    IVec3D ijk(i, j, 0);
    const bool isClassIII = (resolution % 2);
    const int effectiveRes = isClassIII ? resolution + 1 : resolution;

    digits.assign(resolution + 1, IVec3D::INVALID_DIGIT);

    bool first = true;
    for (int r = effectiveRes; r >= 1; r--) {
        IVec3D lastIJK = ijk;
        IVec3D lastCenter;
        if (r % 2) {
            ijk.upAp7();
            lastCenter = ijk;
            lastCenter.downAp7();
        } else {
            ijk.upAp7r();
            lastCenter = ijk;
            lastCenter.downAp7r();
        }
        if (first && isClassIII) {
            first = false;
            continue;
        }
        IVec3D diff = lastIJK.diffVec(lastCenter);
        digits[r] = diff.unitIjkPlusToDigit();
    }
    return ijk;
}

std::string encode_bijective(int quadNum, long long i, long long j, int resolution) {
    if (resolution == 0) {
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << quadNum;
        return oss.str();
    }

    std::vector<IVec3D::Direction> digits;
    const IVec3D::Direction seed =
        z7_walk_up(i, j, resolution, digits).unitIjkPlusToDigit();
    if (seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error(
            "Z7 encode: coordinate does not lie in the given quad");
    }

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << (quadNum + topo().n_quads() * seed);
    std::string out = oss.str();
    for (int r = 1; r <= resolution; r++) {
        out += std::to_string((int) digits[r]);
    }
    return out;
}

void decode_bijective(const std::string& index, int resolution,
                      int& quadNum, long long& i, long long& j) {
    if (index.length() < 2) {
        throw std::runtime_error("Z7 index too short");
    }
    const int lead = std::stoi(index.substr(0, 2));
    const int n_quads = topo().n_quads();
    quadNum = lead % n_quads;
    const int seed = lead / n_quads;
    if (lead < 0 || seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error("Invalid base cell number");
    }

    std::string z7str = index.substr(2);
    int res = (int) z7str.length();
    if (res == 0) {
        // Resolution 0 is one cell per vertex of the solid, so the seed is all
        // that says which of the base cells meeting at the quad's corner the
        // index names. encode_bijective() writes the quad's own, and stripping
        // a digit off a deeper index leaves whichever the ancestry came from.
        quadNum = seed_base_cell(quadNum, seed);
        i = 0;
        j = 0;
        return;
    }
    if (res % 2) {
        z7str += "0";
        res++;
    }

    IVec3D ijk = z7_seed_coord(seed);
    for (int r = 0; r < res; r++) {
        if ((r + 1) % 2) {
            ijk.downAp7();
        } else {
            ijk.downAp7r();
        }
        ijk.neighbor((IVec3D::Direction) (z7str.c_str()[r] - '0'));
    }

    IVec2D ij(ijk);
    i = ij.i();
    j = ij.j();
}

std::string canonical_form(const std::string& z7_index, int max_iterations) {
    // Handle resolution 0 (just base cell)
    if (z7_index.length() <= 2) {
        return z7_index;
    }
    
    std::unordered_set<std::string> seen;
    std::vector<std::string> orbit;
    std::string current = z7_index;
    
    // Iterate until we find a cycle or fixed point
    for (int iter = 0; iter < max_iterations; ++iter) {
        if (seen.count(current)) {
            // Found a cycle - return the lexicographically smallest in the orbit
            return *std::min_element(orbit.begin(), orbit.end());
        }
        
        seen.insert(current);
        orbit.push_back(current);
        
        // Decode and re-encode to get next in sequence
        int quadNum;
        long long i, j;
        int res = current.length() - 2;
        
        decode_bijective(current, res, quadNum, i, j);
        std::string next = encode_bijective(quadNum, i, j, res);
        
        // Check for fixed point
        if (next == current) {
            return current;
        }
        
        current = next;
    }
    
    // If we hit max iterations without finding a cycle, 
    // return the smallest we've seen so far
    return *std::min_element(orbit.begin(), orbit.end());
}

} // namespace z7
} // namespace hexify
