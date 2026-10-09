// index_z7.cpp
// Z7 Hierarchical Index Implementation for Aperture-7 Hexagonal Grids
//
// Encodes (quad, i, j) coordinates as a hierarchical digit label.
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
//   reference implementation, DgZ7System its 64-bit integer form
// - Sahr (2019) "Central Place Indexing", Cartographica 54(1): 16-29, for
//   generalized balanced ternary digit arithmetic
// - H3 coordijk.c and algos.c (Apache 2.0) for aperture-7 coordinate math and
//   neighbours by digit arithmetic
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "index_z7.h"
#include "polyhedron.h"
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <algorithm>

namespace hexify {
namespace z7 {

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
// resolution 0, recording the child ordinal of each level in digit[1..res].
// The return value is the coordinate the walk arrives at: the origin when every
// ancestor of the cell lies in the same quad, and a neighbouring lattice point
// when the quad's rhombic boundary cuts through one of them.
static IVec3D z7_walk_up(long long i, long long j, int resolution, Digit* digit) {
    IVec3D ijk(i, j, 0);
    const bool isClassIII = (resolution % 2);
    const int effectiveRes = isClassIII ? resolution + 1 : resolution;

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
        if (r > resolution) continue;
        digit[r] = lastIJK.diffVec(lastCenter).unitIjkPlusToDigit();
    }
    return ijk;
}

// The lattice point a label's digits name, walked down from `ijk` (padded to
// the Class I substrate at odd resolutions).
static IVec3D z7_walk_down(IVec3D ijk, const Label& label) {
    const int res = label.res;
    for (int r = 1; r <= res; r++) {
        if (r % 2) {
            ijk.downAp7();
        } else {
            ijk.downAp7r();
        }
        ijk.neighbor(label.digit[r]);
    }
    if (res % 2) ijk.downAp7r();
    return ijk;
}

// ============================================================================
// Digit rotations
// ============================================================================

// The first nonzero digit, or the centre digit when there is none.
static inline Digit leading_digit(const Label& label) {
    for (int r = 1; r <= label.res; r++) {
        if (label.digit[r] != IVec3D::CENTER_DIGIT) return label.digit[r];
    }
    return IVec3D::CENTER_DIGIT;
}

static inline void rotate_ccw(Label& label) {
    for (int r = 1; r <= label.res; r++) label.digit[r] = IVec3D::rotate60ccw(label.digit[r]);
}

static inline void rotate_cw(Label& label) {
    for (int r = 1; r <= label.res; r++) label.digit[r] = IVec3D::rotate60cw(label.digit[r]);
}

// A 60-degree counter-clockwise turn about a pentagon whose direction `skip`
// is deleted: a label turned onto that direction turns once more.
static inline void rotate_pent_ccw(Label& label, Digit skip) {
    rotate_ccw(label);
    if (leading_digit(label) == skip) rotate_ccw(label);
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

bool igeo7_labels() { return topo().solid == Solid::Icosahedron; }

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

static void name_igeo7(int quadNum, int seed, Label& label) {
    const SolidTopology& t = topo();
    const int baseCell = seed_base_cell(quadNum, seed);
    const Digit skip = skipped_digit(baseCell);

    if (baseCell != quadNum && t.is_pole(baseCell)) {
        const int turns = northern_base_cell(baseCell) ? quadNum - 1
                                                       : t.n_diamonds() - quadNum;
        for (int k = 0; k < turns; k++) {
            IVec3D::rotateDigitVecCCW(label.digit, label.res, skip);
        }
    }

    label.lead = baseCell;
    if (leading_digit(label) == skip) rotate_ccw(label);
}

static void decode_igeo7(const Label& label, int& quadNum, long long& i, long long& j) {
    const int baseCell = label.lead;
    const int res = label.res;
    quadNum = baseCell;
    i = 0;
    j = 0;
    if (res == 0) return;

    // The quad edge in Class I substrate units.
    long long unit = 1;
    for (int r = 0; r < (res + 1) / 2; r++) unit *= 7;

    IVec2D ij(z7_walk_down(IVec3D(0, 0, 0), label));
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

static void decode_quad_seed(const Label& label, int& quadNum, long long& i, long long& j) {
    const int n_quads = topo().n_quads();
    quadNum = label.lead % n_quads;
    const int seed = label.lead / n_quads;
    if (seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error("Invalid base cell number");
    }
    if (label.res == 0) {
        // Resolution 0 is one cell per vertex of the solid, so the seed is all
        // that says which of the base cells meeting at the quad's corner the
        // index names.
        quadNum = seed_base_cell(quadNum, seed);
        i = 0;
        j = 0;
        return;
    }
    IVec2D ij(z7_walk_down(z7_seed_coord(seed), label));
    i = ij.i();
    j = ij.j();
}

// ============================================================================
// Labels
// ============================================================================

Label label_of(int quadNum, long long i, long long j, int resolution) {
    if (resolution < 0 || resolution > kMaxResolution) {
        throw std::runtime_error("Z7: resolution out of range");
    }
    Label label;
    label.res = resolution;
    const Digit seed = z7_walk_up(i, j, resolution, label.digit).unitIjkPlusToDigit();
    if (seed >= IVec3D::NUM_DIGITS) {
        throw std::runtime_error(
            "Z7 encode: coordinate does not lie in the given quad");
    }
    if (igeo7_labels()) {
        name_igeo7(quadNum, seed, label);
    } else {
        label.lead = quadNum + topo().n_quads() * seed;
    }
    return label;
}

void cell_of(const Label& label, int& quadNum, long long& i, long long& j) {
    if (igeo7_labels()) {
        decode_igeo7(label, quadNum, i, j);
    } else {
        decode_quad_seed(label, quadNum, i, j);
    }
}

std::string to_string(const Label& label) {
    std::string out(2 + label.res, '0');
    out[0] = (char) ('0' + label.lead / 10);
    out[1] = (char) ('0' + label.lead % 10);
    for (int r = 1; r <= label.res; r++) out[1 + r] = (char) ('0' + (int) label.digit[r]);
    return out;
}

Label from_string(const std::string& index) {
    if (index.length() < 2) {
        throw std::runtime_error("Z7 index too short");
    }
    const int res = (int) index.length() - 2;
    if (res > kMaxResolution) {
        throw std::runtime_error("Z7 index longer than the deepest resolution");
    }
    Label label;
    if (index[0] < '0' || index[0] > '9' || index[1] < '0' || index[1] > '9') {
        throw std::runtime_error("Invalid base cell number");
    }
    label.lead = 10 * (index[0] - '0') + (index[1] - '0');
    if (igeo7_labels() && label.lead >= topo().n_quads()) {
        throw std::runtime_error("Invalid base cell number");
    }
    label.res = res;
    for (int r = 1; r <= res; r++) {
        const int d = index[1 + r] - '0';
        if (d < 0 || d > 6) throw std::runtime_error("Z7 digits must be 0-6");
        label.digit[r] = (Digit) d;
    }
    return label;
}

// ============================================================================
// Packed and monotonic forms (IGEO7)
// ============================================================================

static inline int packed_shift(int r) { return 3 * (kMaxPackedRes - r); }

uint64_t to_packed(const Label& label) {
    if (!igeo7_labels()) {
        throw std::runtime_error("Z7: the packed form is defined for IGEO7 labels only");
    }
    if (label.res > kMaxPackedRes) {
        throw std::runtime_error("Z7: the packed form holds resolutions 0 to 20");
    }
    uint64_t z = static_cast<uint64_t>(label.lead) << 60;
    for (int r = 1; r <= kMaxPackedRes; r++) {
        const uint64_t d = (r <= label.res) ? static_cast<uint64_t>(label.digit[r])
                                            : static_cast<uint64_t>(IVec3D::INVALID_DIGIT);
        z |= d << packed_shift(r);
    }
    return z;
}

bool from_packed(uint64_t z, Label& out) {
    out.lead = static_cast<int>(z >> 60);
    if (out.lead >= topo().n_quads()) return false;
    out.res = kMaxPackedRes;
    for (int r = 1; r <= kMaxPackedRes; r++) {
        const int d = static_cast<int>((z >> packed_shift(r)) & 7u);
        if (d == IVec3D::INVALID_DIGIT) {
            if (out.res == kMaxPackedRes) out.res = r - 1;
        } else if (out.res != kMaxPackedRes) {
            return false;
        } else {
            out.digit[r] = (Digit) d;
        }
    }
    return true;
}

uint64_t to_monotonic(const Label& label) {
    if (label.res > kMaxMonotonicRes) {
        throw std::runtime_error("Z7: monotonic IDs pass 2^63 beyond resolution 21");
    }
    uint64_t m = static_cast<uint64_t>(label.lead);
    for (int r = 1; r <= label.res; r++) m = 7 * m + static_cast<uint64_t>(label.digit[r]);
    return m;
}

bool from_monotonic(uint64_t m, int resolution, Label& out) {
    if (resolution < 0 || resolution > kMaxMonotonicRes) return false;
    out.res = resolution;
    for (int r = resolution; r >= 1; r--) {
        out.digit[r] = (Digit) (m % 7);
        m /= 7;
    }
    if (m >= static_cast<uint64_t>(topo().n_quads())) return false;
    out.lead = static_cast<int>(m);
    return true;
}

// ============================================================================
// Neighbours by digit arithmetic (generalized balanced ternary)
// ============================================================================
// A cell at resolution r is the lattice point sum over levels of its digits'
// unit vectors, each scaled down through the levels below it. Adding a unit
// vector at level r changes digit r and may carry a unit vector into level
// r - 1: a + d = new digit + carry scaled down one level. The sum and carry
// tables follow from the level's scaling, downAp7 at odd levels and downAp7r
// at even ones, so they are built from the same operations the walks use.

struct DigitSumTables {
    Digit sum[2][7][7];    // [level is odd][digit][direction]
    Digit carry[2][7][7];
};

static DigitSumTables make_digit_sum_tables() {
    DigitSumTables t;
    for (int odd = 0; odd < 2; odd++) {
        for (int a = 0; a < 7; a++) {
            for (int d = 0; d < 7; d++) {
                IVec3D target = z7_seed_coord(a);
                target.neighbor((Digit) d);
                bool found = false;
                for (int c = 0; c < 7 && !found; c++) {
                    for (int n = 0; n < 7 && !found; n++) {
                        IVec3D v = z7_seed_coord(c);
                        if (odd) v.downAp7(); else v.downAp7r();
                        v.neighbor((Digit) n);
                        if (v == target) {
                            t.sum[odd][a][d] = (Digit) n;
                            t.carry[odd][a][d] = (Digit) c;
                            found = true;
                        }
                    }
                }
                if (!found) throw std::runtime_error("Z7: digit sum table has no entry");
            }
        }
    }
    return t;
}

static const DigitSumTables& digit_sums() {
    static const DigitSumTables t = make_digit_sum_tables();
    return t;
}

// IGEO7 base-cell adjacency in the base cells' label frames, indexed by base
// cell and direction 1-6. A digit sum that carries direction c out of
// resolution 1 of base cell b lands in base cell kBaseNeighbor[b][c], and the
// digits, still in b's frame, turn kBaseTurns[b][c] times by 60 degrees
// counter-clockwise about the new pentagon into its frame, each turn passing
// over its deleted direction. Base cells 1-5 form the northern ring and 6-10
// the southern one; a ring cell reaches its pole in its deleted direction and
// the direction beside it, and the pole's turns go round with the ring. The
// tables play the role of H3's baseCellNeighbors and
// baseCellNeighbor60CCWRots for IGEO7's twelve pentagons.
static const int kBaseNeighbor[12][7] = {
    {-1,  5,  4,  4,  2,  1,  3},
    {-1,  5,  0,  0,  6, 10,  2},
    {-1,  1,  0,  0,  7,  6,  3},
    {-1,  2,  0,  0,  8,  7,  4},
    {-1,  3,  0,  0,  9,  8,  5},
    {-1,  4,  0,  0, 10,  9,  1},
    {-1, 10,  2,  1, 11, 11,  7},
    {-1,  6,  3,  2, 11, 11,  8},
    {-1,  7,  4,  3, 11, 11,  9},
    {-1,  8,  5,  4, 11, 11, 10},
    {-1,  9,  1,  5, 11, 11,  6},
    {-1,  9,  6, 10,  8,  8,  7}};

static const int kBaseTurns[12][7] = {
    {0, 1, 3, 2, 5, 0, 4},
    {0, 0, 0, 4, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 0},
    {0, 0, 2, 1, 0, 0, 0},
    {0, 0, 3, 2, 0, 0, 0},
    {0, 0, 4, 3, 0, 0, 0},
    {0, 0, 0, 0, 4, 0, 0},
    {0, 0, 0, 0, 3, 4, 0},
    {0, 0, 0, 0, 2, 3, 0},
    {0, 0, 0, 0, 1, 2, 0},
    {0, 0, 0, 0, 0, 1, 0},
    {0, 4, 0, 5, 2, 3, 1}};

bool neighbor(const Label& in, Digit dir, Label& out) {
    const DigitSumTables& t = digit_sums();
    out = in;
    Digit carry = dir;
    for (int r = in.res; r >= 1 && carry != IVec3D::CENTER_DIGIT; r--) {
        const int odd = r & 1;
        const Digit a = out.digit[r];
        out.digit[r] = t.sum[odd][a][carry];
        carry = t.carry[odd][a][carry];
    }

    if (carry != IVec3D::CENTER_DIGIT) {
        out.lead = kBaseNeighbor[in.lead][carry];
        const Digit skip = skipped_digit(out.lead);
        for (int k = 0; k < kBaseTurns[in.lead][carry]; k++) rotate_pent_ccw(out, skip);
        if (leading_digit(out) == skip) rotate_ccw(out);
        return true;
    }

    // Inside the base cell, a step into the deleted subsequence crosses the
    // pentagon's seam: the cells on its two sides are a 60-degree turn apart.
    const Digit skip = skipped_digit(in.lead);
    if (leading_digit(out) != skip) return true;
    const Digit from = leading_digit(in);
    if (from == IVec3D::CENTER_DIGIT) return false;
    if (from == IVec3D::rotate60ccw(skip)) {
        rotate_cw(out);
    } else if (from == IVec3D::rotate60cw(skip)) {
        rotate_ccw(out);
    } else {
        throw std::runtime_error("Z7: a neighbour step reached the deleted "
                                 "subsequence from a direction not beside it");
    }
    return true;
}

// ============================================================================

std::string encode(int quadNum, long long i, long long j, int resolution) {
    return to_string(label_of(quadNum, i, j, resolution));
}

void decode(const std::string& index, int& quadNum, long long& i, long long& j) {
    cell_of(from_string(index), quadNum, i, j);
}

bool in_deleted_subsequence(const Label& label) {
    if (!igeo7_labels()) return false;
    const Digit lead = leading_digit(label);
    return lead != IVec3D::CENTER_DIGIT && lead == skipped_digit(label.lead);
}

bool in_deleted_subsequence(const std::string& index) {
    if (!igeo7_labels() || index.length() < 2) return false;
    return in_deleted_subsequence(from_string(index));
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
