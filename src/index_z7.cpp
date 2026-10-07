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
// 3. Name the base cell the walk arrives in, and turn the digits into that
//    base cell's frame (on the icosahedron, IGEO7's Z7 index)
//
// References:
// - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
// - Kmoch, Sahr, Chan, Uuemaa (2025) "IGEO7: A new hierarchically indexed
//   hexagonal equal-area discrete global grid system", AGILE: GIScience
//   Series 6, doi:10.5194/agile-giss-6-32-2025; DGGRID's DgZ7StringRF is its
//   reference implementation
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

typedef IVec3D::Direction Digit;

static std::string two_digits(int value) {
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << value;
    return oss.str();
}

// The level-0 coordinate a walk arrives at. A cell whose whole ancestry lies
// inside its quad arrives at the origin. The quad is a rhombus while the
// aperture-7 parents are hexagons, so the quad boundary cuts through the
// parents of the cells along it; those arrive at one of the six neighbours of
// the origin instead, a unit digit.
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
// other three corners, whose vertex is the base cell the cell descends from.
// A vertex quad has no other corners.
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
                         std::vector<Digit>& digits) {
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

// The lattice point a digit string names, walked down from `start` (padded to
// the Class I substrate at odd resolutions).
static IVec3D z7_walk_down(IVec3D ijk, std::string z7str) {
    int res = (int) z7str.length();
    if (res % 2) {
        z7str += "0";
        res++;
    }
    for (int r = 0; r < res; r++) {
        if ((r + 1) % 2) {
            ijk.downAp7();
        } else {
            ijk.downAp7r();
        }
        ijk.neighbor((Digit) (z7str[r] - '0'));
    }
    return ijk;
}

// ============================================================================
// IGEO7 labels (icosahedron)
// ============================================================================
// Every resolution-0 cell of the icosahedron is a pentagon centred on a
// vertex, and IGEO7 names a cell by the base cell it descends from and the
// digits read in that base cell's frame. Three steps turn a walk in a quad's
// frame into that label, as DGGRID's DgZ7StringRF::quantify does:
//
// - a cell whose walk leaves the quad is named under the base cell at the
//   corner it arrives at (seed_base_cell());
// - under a polar base cell, the digits of a cell from the k-th quad around
//   the pole are turned k times by 60 degrees into the pole's frame;
// - a pentagon lacks one of the six directions, so the subsequence whose first
//   nonzero digit is the missing one is deleted: such a cell's digits are
//   turned once more, onto the next direction, which is free in that frame.
//   The missing direction is digit 2 under the five northern base cells and
//   the north pole, and digit 5 under the others.
//
// Decoding unfolds the lattice point a label names back into the quad that
// owns it, as DgZ7StringRF::invQuantify does. The quads around the north pole
// are 1..n_ring and those around the south pole n_ring+1..2*n_ring.

static int ring_size() { return topo().n_diamonds() / 2; }

static bool northern_base_cell(int baseCell) { return baseCell <= ring_size(); }

static Digit skipped_digit(int baseCell) {
    return (Digit) (northern_base_cell(baseCell) ? IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE1
                                                 : IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE2);
}

// The diamond quad whose corner `corner` is vertex `vertex`.
static int quad_with_corner(int corner, int vertex) {
    const SolidTopology& t = topo();
    for (int q = 1; q <= t.n_diamonds(); q++) {
        if (t.corner[q][corner] == vertex) return q;
    }
    throw std::runtime_error("Z7: no quad has the requested corner");
}

static std::string encode_igeo7(int quadNum, int seed, std::vector<Digit>& digits,
                                int resolution) {
    const SolidTopology& t = topo();
    const int baseCell = seed_base_cell(quadNum, seed);
    const Digit skip = skipped_digit(baseCell);

    if (baseCell != quadNum && t.is_pole(baseCell)) {
        const int turns = northern_base_cell(baseCell) ? quadNum - 1
                                                       : t.n_diamonds() - quadNum;
        for (int k = 0; k < turns; k++) {
            IVec3D::rotateDigitVecCCW(digits.data(), resolution, skip);
        }
    }

    std::string out = two_digits(baseCell);
    bool seenNonZero = false;
    bool skipRotate = false;
    for (int r = 1; r <= resolution; r++) {
        Digit d = digits[r];
        if (!seenNonZero && d != IVec3D::CENTER_DIGIT) {
            seenNonZero = true;
            skipRotate = (d == skip);
        }
        if (skipRotate) d = IVec3D::rotate60ccw(d);
        out += (char) ('0' + (int) d);
    }
    return out;
}

static void decode_igeo7(int baseCell, const std::string& z7str,
                         int& quadNum, long long& i, long long& j) {
    const int res = (int) z7str.length();
    quadNum = baseCell;
    i = 0;
    j = 0;
    if (res == 0) return;

    // The quad edge in Class I substrate units.
    long long unit = 1;
    for (int r = 0; r < (res + 1) / 2; r++) unit *= 7;

    IVec2D ij(z7_walk_down(IVec3D(0, 0, 0), z7str));
    if (ij.i() == 0 && ij.j() == 0) return;

    const int n = ring_size();
    const SolidTopology& t = topo();
    const bool negI = ij.i() < 0;
    const bool negJ = ij.j() < 0;
    const long long origI = ij.i();

    if (baseCell == 0) {
        if (!negI) {
            if (!negJ) {
                if (ij.i() > ij.j()) {
                    quadNum = 2;
                    ij.setI(ij.j());
                    ij.setJ(unit - (origI - ij.j()));
                } else {
                    quadNum = 3;
                    ij.setI(ij.j() - ij.i());
                    ij.setJ(unit - origI);
                }
            } else {
                quadNum = 1;
                ij.setJ(ij.j() + unit);
            }
        } else {
            if (!negJ) {
                if (ij.j() == 0) {
                    quadNum = 4;
                    ij.setJ(unit + ij.i());
                    ij.setI(0);
                } else {
                    quadNum = 3;
                    ij.setI(-ij.i());
                    ij.setJ(unit - ij.j());
                }
            } else {
                if (ij.i() < ij.j()) {
                    quadNum = 4;
                    ij.setI(-ij.j());
                    ij.setJ(unit - (-origI + ij.j()));
                } else {
                    quadNum = 5;
                    ij.setI(origI - ij.j());
                    ij.setJ(unit + origI);
                }
            }
        }
    } else if (baseCell == t.south_pole()) {
        if (!negI) {
            if (!negJ) {
                if (ij.i() == 0) {
                    quadNum = n + 1;
                    ij.setI(unit - ij.j());
                    ij.setJ(0);
                } else if (ij.j() == 0) {
                    quadNum = n + 3;
                    ij.setI(unit - ij.i());
                    ij.setJ(0);
                } else if (ij.j() > ij.i()) {
                    quadNum = n + 1;
                    ij.setI(unit - (ij.j() - ij.i()));
                    ij.setJ(origI);
                } else {
                    quadNum = n + 2;
                    ij.setI(unit - ij.j());
                    ij.setJ(origI - ij.j());
                }
            } else {
                quadNum = n + 3;
                ij.setI(unit - ij.i());
                ij.setJ(-ij.j());
            }
        } else {
            if (negJ) {
                if (ij.i() > ij.j()) {
                    quadNum = n + 3;
                    ij.setI(unit - (-ij.j() + ij.i()));
                    ij.setJ(-origI);
                } else {
                    quadNum = n + 4;
                    ij.setI(unit + ij.j());
                    ij.setJ(-origI + ij.j());
                }
            } else {
                quadNum = n + 5;
                ij.setI(unit + ij.i());
            }
        }
    } else if (northern_base_cell(baseCell)) {
        if (negJ) {
            ij.setJ(ij.j() + unit);
            if (negI) {
                ij.setI(ij.i() + unit);
                quadNum = quad_with_corner(kCornerFar, baseCell);
            } else {
                quadNum = quad_with_corner(kCornerJ, baseCell);
            }
        } else if (negI) {
            // the deleted direction 2: digit 3 was turned onto it
            IVec3D ijk(ij.i(), ij.j(), 0);
            ijk.ijkRotate60cw();
            ij = IVec2D(ijk);
        }
    } else {
        if (negI) {
            ij.setI(ij.i() + unit);
            if (negJ) {
                ij.setJ(ij.j() + unit);
                quadNum = quad_with_corner(kCornerFar, baseCell);
            } else {
                quadNum = quad_with_corner(kCornerI, baseCell);
            }
        } else if (negJ) {
            ij.setI(ij.j() + unit);
            ij.setJ(ij.j() + unit - origI);
            quadNum = quad_with_corner(kCornerFar, baseCell);
        }
    }

    i = ij.i();
    j = ij.j();
}

// ============================================================================
// Quad-and-seed labels (other solids)
// ============================================================================
// IGEO7 is defined on the icosahedron, whose vertex cells are pentagons; the
// octahedron's are squares. There the quad stays fixed and the leading field
// carries the arrival point as quad + n * seed, n the solid's number of quads:
// the plain quad whenever the cell's ancestry stays inside it.

static std::string encode_quad_seed(int quadNum, int seed,
                                    const std::vector<Digit>& digits, int resolution) {
    std::string out = two_digits(quadNum + topo().n_quads() * seed);
    for (int r = 1; r <= resolution; r++) {
        out += (char) ('0' + (int) digits[r]);
    }
    return out;
}

static void decode_quad_seed(int lead, const std::string& z7str,
                             int& quadNum, long long& i, long long& j) {
    const int n_quads = topo().n_quads();
    quadNum = lead % n_quads;
    const int seed = lead / n_quads;
    if (seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error("Invalid base cell number");
    }
    if (z7str.empty()) {
        // Resolution 0 is one cell per vertex of the solid, so the seed is all
        // that says which of the base cells meeting at the quad's corner the
        // index names.
        quadNum = seed_base_cell(quadNum, seed);
        i = 0;
        j = 0;
        return;
    }
    IVec2D ij(z7_walk_down(z7_seed_coord(seed), z7str));
    i = ij.i();
    j = ij.j();
}

// ============================================================================

static bool igeo7_labels() { return topo().solid == Solid::Icosahedron; }

std::string encode(int quadNum, long long i, long long j, int resolution) {
    std::vector<Digit> digits;
    const Digit seed = z7_walk_up(i, j, resolution, digits).unitIjkPlusToDigit();
    if (seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error(
            "Z7 encode: coordinate does not lie in the given quad");
    }
    return igeo7_labels() ? encode_igeo7(quadNum, seed, digits, resolution)
                          : encode_quad_seed(quadNum, seed, digits, resolution);
}

void decode(const std::string& index, int& quadNum, long long& i, long long& j) {
    if (index.length() < 2) {
        throw std::runtime_error("Z7 index too short");
    }
    const int lead = std::stoi(index.substr(0, 2));
    if (lead < 0 || (igeo7_labels() && lead >= topo().n_quads())) {
        throw std::runtime_error("Invalid base cell number");
    }
    const std::string z7str = index.substr(2);
    if (igeo7_labels()) {
        decode_igeo7(lead, z7str, quadNum, i, j);
    } else {
        decode_quad_seed(lead, z7str, quadNum, i, j);
    }
}

bool in_deleted_subsequence(const std::string& index) {
    if (!igeo7_labels() || index.length() < 2) return false;
    const Digit skip = skipped_digit(std::stoi(index.substr(0, 2)));
    for (size_t r = 2; r < index.length(); r++) {
        if (index[r] != '0') return index[r] - '0' == (int) skip;
    }
    return false;
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

        decode(current, quadNum, i, j);
        std::string next = encode(quadNum, i, j, res);

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
