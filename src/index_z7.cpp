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
// 3. Handle pentagon base cells (0, 11) with rotation adjustments
// 4. Encode digits 0-6 representing the 7 child positions
//
// References:
// - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
// - H3 coordijk.c (Apache 2.0) for aperture-7 coordinate math
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "index_z7.h"
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <vector>
#include <algorithm>

namespace hexify {
namespace z7 {

const int adjacentBaseCellTable[12][4] = {
    { 0, 0, 0, 0 },
    { 1, 6, 2, 0 },
    { 2, 7, 3, 0 },
    { 3, 8, 4, 0 },
    { 4, 9, 5, 0 },
    { 5, 10, 1, 0 },
    { 6, 11, 7, 2 },
    { 7, 11, 8, 3 },
    { 8, 11, 9, 4 },
    { 9, 11, 10, 5 },
    { 10, 11, 6, 1 },
    { 11, 11, 0, 0 }
};

const int inverseAdjacentBaseCellTable[12][2] = {
    { 0,  0 },
    { 5, 10 },
    { 1,  6 },
    { 2,  7 },
    { 3,  8 },
    { 4,  9 },
    { 10, 1 },
    { 6,  2 },
    { 7,  3 },
    { 8,  4 },
    { 9,  5 },
    { 11, 11 }
};

std::string encode(int quadNum, long long i, long long j, int resolution) {
    // Format base cell as 2 digits
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << quadNum;
    std::string bcstr = oss.str();
    
    // Resolution 0 is just the base cell
    if (resolution == 0) {
        return bcstr;
    }
    
    IVec3D ijk(i, j, 0);
    int baseCell = quadNum;
    IVec3D baseCellIjk = ijk;
    int res = resolution;
    
    bool isClassIII = (res % 2);
    int effectiveRes = (isClassIII) ? res + 1 : res;
    
    // Allocate (res + 1) elements for digit storage
    // Index 0 unused; indices 1..res store the hierarchical path digits
    IVec3D::Direction* digits = 
        (IVec3D::Direction*) malloc((res + 1) * sizeof(IVec3D::Direction));
    for (int r = 0; r < res + 1; r++) {
        digits[r] = IVec3D::INVALID_DIGIT;
    }
    
    bool first = true;
    for (int r = effectiveRes; r >= 0; r--) {
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
        
        if (r == 1) {
            baseCellIjk = ijk;
        }
        
        if (first && isClassIII) {
            first = false;
            continue;
        }
        
        IVec3D diff = lastIJK.diffVec(lastCenter);
        digits[r] = diff.unitIjkPlusToDigit();
    }
    
    int quadOriginBaseCell = baseCell;
    
    // Apply adjacency transformations
    if (baseCellIjk.i() == 1) {
        if (baseCellIjk.j() == 0) {
            baseCell = adjacentBaseCellTable[baseCell][1];
        } else {
            baseCell = adjacentBaseCellTable[baseCell][2];
        }
    } else if (baseCellIjk.j() == 1) {
        baseCell = adjacentBaseCellTable[baseCell][3];
    }
    
    // Handle the single-cell quads 0 and 11
    if (baseCell != quadOriginBaseCell) {
        if (baseCell == 0) {
            // must be quad 1 - 5
            // rotate once for each quad past 1
            for (int q = 1; q < quadOriginBaseCell; q++) {
                IVec3D::rotateDigitVecCCW(digits, res, 
                    (IVec3D::Direction)IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE1);
            }
        } else if (baseCell == 11) {
            // must be quad 6 - 10
            // rotate once for each quad less than 10
            int numRots = 10 - quadOriginBaseCell;
            for (int q = 0; q < numRots; q++) {
                IVec3D::rotateDigitVecCCW(digits, res,
                    (IVec3D::Direction)IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE2);
            }
        }
    }
    
    // Format the base cell for output
    oss.str("");
    oss << std::setfill('0') << std::setw(2) << baseCell;
    std::string addstr = oss.str();
    
    // Pentagon digit skip handling for base cells near poles
    // Pentagons have only 5 neighbors (not 6), so one digit direction is invalid.
    // Type 1 (north pole region, cells 0-5): skip digit 2 (J_AXES_DIGIT)
    // Type 2 (south pole region, cells 6-11): skip digit 5 (IJ_AXES_DIGIT)
    IVec3D::Direction skipDigit = ((baseCell < 6) ?
        (IVec3D::Direction)IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE1 :
        (IVec3D::Direction)IVec3D::PENTAGON_SKIPPED_DIGIT_TYPE2);

    // Rotation adjustment for pentagon continuity
    // When the first non-zero digit equals the skipped digit, rotate all
    // subsequent digits to maintain consistent indexing across the pentagon gap.
    bool skipRotate = false;
    bool firstNonZero = false;

    for (int r = 1; r < res + 1; r++) {
        IVec3D::Direction d = digits[r];

        // Check if this is the first non-zero digit
        if (!firstNonZero && d != IVec3D::CENTER_DIGIT) {
            firstNonZero = true;
            if (d == skipDigit)
                skipRotate = true;
        }

        // Apply rotation to all digits when skip condition is triggered
        if (skipRotate) {
            d = IVec3D::rotate60ccw(d);
        }
        
        addstr += std::to_string((int)d);
    }
    
    free(digits);
    
    return addstr;
}

void decode(const std::string& z7_index, int resolution,
            int& quadNum, long long& i, long long& j) {
    
    if (z7_index.length() < 2) {
        throw std::runtime_error("Z7 index too short");
    }
    
    std::string bcStr = z7_index.substr(0, 2);
    if (bcStr[0] == '0') {
        bcStr = bcStr.substr(1, 1);
    }
    int bcNum = std::stoi(bcStr);
    
    if (bcNum < 0 || bcNum > 11) {
        throw std::runtime_error("Invalid base cell number");
    }
    
    std::string z7str = z7_index.substr(2);
    int res = (int) z7str.length();
    
    // Resolution 0 is just the base cell
    if (res == 0) {
        quadNum = bcNum;
        i = 0;
        j = 0;
        return;
    }
    
    if (res % 2) {
        z7str += "0";
        res++;
    }
    
    IVec3D ijk(0, 0, 0);
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
    quadNum = bcNum;
    
    if (i == 0 && j == 0) {
        return;
    }
    
    int numClassI = (resolution + 1) / 2;
    unsigned long long int unitScaleClassIres = 1;
    for (int r = 0; r < numClassI; r++) {
        unitScaleClassIres *= 7;
    }
    
    bool negI = (i < 0);
    bool negJ = (j < 0);
    
    long long origI = i;
    
    // Apply adjacency transformations
    if (bcNum == 0) {
        if (!negI) {
            if (!negJ) {
                if (i > j) {
                    quadNum = 2;
                    i = j;
                    j = unitScaleClassIres - (origI - j);
                } else {
                    quadNum = 3;
                    i = j - i;
                    j = unitScaleClassIres - origI;
                }
            } else {
                quadNum = 1;
                j = j + unitScaleClassIres;
            }
        } else {
            if (!negJ) {
                if (j == 0) {
                    quadNum = 4;
                    j = unitScaleClassIres + i;
                    i = 0;
                } else {
                    quadNum = 3;
                    i = -i;
                    j = unitScaleClassIres - j;
                }
            } else {
                if (i < j) {
                    quadNum = 4;
                    i = -j;
                    j = unitScaleClassIres - (-origI + j);
                } else {
                    quadNum = 5;
                    i = origI - j;
                    j = unitScaleClassIres + origI;
                }
            }
        }
    }
    else if (bcNum == 11) {
        if (!negI) {
            if (!negJ) {
                if (i == 0) {
                    quadNum = 6;
                    i = unitScaleClassIres - j;
                    j = 0;
                } else if (j == 0) {
                    quadNum = 8;
                    i = unitScaleClassIres - i;
                    j = 0;
                } else if (j > i) {
                    quadNum = 6;
                    i = unitScaleClassIres - (j - i);
                    j = origI;
                } else {
                    quadNum = 7;
                    i = unitScaleClassIres - j;
                    j = origI - j;
                }
            } else {
                quadNum = 8;
                i = unitScaleClassIres - i;
                j = -j;
            }
        } else {
            if (negJ) {
                if (i > j) {
                    quadNum = 8;
                    i = unitScaleClassIres - (-j + i);
                    j = -origI;
                } else {
                    quadNum = 9;
                    i = unitScaleClassIres + j;
                    j = -origI + j;
                }
            } else {
                quadNum = 10;
                i = unitScaleClassIres + i;
            }
        }
    }
    else if (bcNum < 6) {
        if (negJ) {
            j = j + unitScaleClassIres;
            if (negI) {
                i = i + unitScaleClassIres;
                quadNum = inverseAdjacentBaseCellTable[bcNum][0];
            } else {
                quadNum = inverseAdjacentBaseCellTable[bcNum][1];
            }
        } else if (negI) {
            IVec3D ijk_temp(i, j);
            ijk_temp.ijkRotate60cw();
            IVec2D ij_temp(ijk_temp);
            i = ij_temp.i();
            j = ij_temp.j();
        }
    }
    else {
        if (negI) {
            i = i + unitScaleClassIres;
            if (negJ) {
                j = j + unitScaleClassIres;
                quadNum = inverseAdjacentBaseCellTable[bcNum][0];
            } else {
                quadNum = inverseAdjacentBaseCellTable[bcNum][1];
            }
        } else if (negJ) {
            i = j + unitScaleClassIres;
            j = j + unitScaleClassIres - origI;
            
            quadNum = inverseAdjacentBaseCellTable[bcNum][0];
        }
    }
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
        
        decode(current, res, quadNum, i, j);
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
