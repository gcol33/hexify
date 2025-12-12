// ijk_coordinates.cpp
// Cube Coordinate Operations for Hexagonal Grids
//
// Implements standard operations on cube coordinates (i,j,k):
// - Normalization: Ensures minimal non-negative representation
// - Rotation: 60° CW/CCW rotations via coordinate cycling
// - Aperture-7 scaling: sqrt(7) factor with ~19.1° rotation
//
// Mathematical basis for aperture-7 transforms:
//   upAp7 (coarsen):   [i', j'] = [3i - j, i + 2j] / 7
//   upAp7r (coarsen):  [i', j'] = [2i + j, 3j - i] / 7
//   downAp7 (refine):  Uses basis vectors {3,0,1}, {1,3,0}, {0,1,3}
//   downAp7r (refine): Uses basis vectors {3,1,0}, {0,3,1}, {1,0,3}
//
// 60° rotation in cube coords: cycle coordinates and negate
//   CCW: (i,j,k) -> (-k,-i,-j) then normalize
//   CW:  (i,j,k) -> (-j,-k,-i) then normalize
//
// References:
// - Red Blob Games "Hexagonal Grids" (cube coordinate operations)
// - H3 coordijk.c (Apache 2.0) for aperture-7 math
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "ijk_coordinates.h"
#include <algorithm>
#include <cmath>

namespace hexify {
namespace z7 {

static const IVec3D UNIT_VECS[] = {
    {0, 0, 0},  // direction 0 (CENTER_DIGIT)
    {0, 0, 1},  // direction 1 (K_AXES_DIGIT)
    {0, 1, 0},  // direction 2 (J_AXES_DIGIT)
    {0, 1, 1},  // direction 3 (I_AXES_DIGIT)
    {1, 0, 0},  // direction 4 (IK_AXES_DIGIT)
    {1, 0, 1},  // direction 5 (IJ_AXES_DIGIT)
    {1, 1, 0}   // direction 6 (JK_AXES_DIGIT)
};

void IVec3D::ijkPlusNormalize() {
    if (i_ < 0) {
        j_ -= i_;
        k_ -= i_;
        i_ = 0;
    }

    if (j_ < 0) {
        i_ -= j_;
        k_ -= j_;
        j_ = 0;
    }

    if (k_ < 0) {
        i_ -= k_;
        j_ -= k_;
        k_ = 0;
    }

    long long min_val = i_;
    if (j_ < min_val) min_val = j_;
    if (k_ < min_val) min_val = k_;
    if (min_val > 0) {
        i_ -= min_val;
        j_ -= min_val;
        k_ -= min_val;
    }
}

IVec3D::Direction IVec3D::unitIjkPlusToDigit() const {
    IVec3D c = *this;
    c.ijkPlusNormalize();

    Direction digit = INVALID_DIGIT;
    for (int i = CENTER_DIGIT; i < NUM_DIGITS; i++) {
        if (c == UNIT_VECS[i]) {
            digit = (Direction) i;
            break;
        }
    }

    return digit;
}

void IVec3D::upAp7() {
    long long i = i_ - k_;
    long long j = j_ - k_;

    i_ = (long long)lround((3 * i - j) / 7.0);
    j_ = (long long)lround((i + 2 * j) / 7.0);
    k_ = 0;
    ijkPlusNormalize();
}

void IVec3D::upAp7r() {
    long long i = i_ - k_;
    long long j = j_ - k_;

    i_ = (long long)lround((2 * i + j) / 7.0);
    j_ = (long long)lround((3 * j - i) / 7.0);
    k_ = 0;
    ijkPlusNormalize();
}

void IVec3D::downAp7() {
    IVec3D iVec = {3, 0, 1};
    IVec3D jVec = {1, 3, 0};
    IVec3D kVec = {0, 1, 3};

    iVec *= i_;
    jVec *= j_;
    kVec *= k_;
    *this = iVec + jVec + kVec;
    
    ijkPlusNormalize();
}

void IVec3D::downAp7r() {
    IVec3D iVec = {3, 1, 0};
    IVec3D jVec = {0, 3, 1};
    IVec3D kVec = {1, 0, 3};

    iVec *= i_;
    jVec *= j_;
    kVec *= k_;

    *this = iVec + jVec + kVec;
    ijkPlusNormalize();
}

void IVec3D::neighbor(Direction digit) {
    if (digit > CENTER_DIGIT && digit < NUM_DIGITS) {
        *this += UNIT_VECS[digit];
        ijkPlusNormalize();
    }
}

void IVec3D::ijkRotate60ccw() {
    IVec3D iVec = {1, 1, 0};
    IVec3D jVec = {0, 1, 1};
    IVec3D kVec = {1, 0, 1};

    iVec *= i_;
    jVec *= j_;
    kVec *= k_;

    *this = iVec + jVec + kVec;
    ijkPlusNormalize();
}

void IVec3D::ijkRotate60cw() {
    IVec3D iVec = {1, 0, 1};
    IVec3D jVec = {1, 1, 0};
    IVec3D kVec = {0, 1, 1};

    iVec *= i_;
    jVec *= j_;
    kVec *= k_;

    *this = iVec + jVec + kVec;
    ijkPlusNormalize();
}

IVec3D::Direction IVec3D::rotate60ccw(Direction digit) {
    switch (digit) {
        case K_AXES_DIGIT:
            return IK_AXES_DIGIT;
        case IK_AXES_DIGIT:
            return I_AXES_DIGIT;
        case I_AXES_DIGIT:
            return IJ_AXES_DIGIT;
        case IJ_AXES_DIGIT:
            return J_AXES_DIGIT;
        case J_AXES_DIGIT:
            return JK_AXES_DIGIT;
        case JK_AXES_DIGIT:
            return K_AXES_DIGIT;
        default:
            return digit;
    }
}

IVec3D::Direction IVec3D::rotate60cw(Direction digit) {
    switch (digit) {
        case K_AXES_DIGIT:
            return JK_AXES_DIGIT;
        case JK_AXES_DIGIT:
            return J_AXES_DIGIT;
        case J_AXES_DIGIT:
            return IJ_AXES_DIGIT;
        case IJ_AXES_DIGIT:
            return I_AXES_DIGIT;
        case I_AXES_DIGIT:
            return IK_AXES_DIGIT;
        case IK_AXES_DIGIT:
            return K_AXES_DIGIT;
        default:
            return digit;
    }
}

void IVec3D::rotateDigitVecCCW(Direction* digits, int maxRes, Direction skipDigit) {
    int foundFirstNonZeroDigit = 0;
    for (int r = 1; r <= maxRes; r++) {
        digits[r] = rotate60ccw(digits[r]);

        if (!foundFirstNonZeroDigit && digits[r] != CENTER_DIGIT) {
            foundFirstNonZeroDigit = 1;

            if (digits[r] == skipDigit)
                rotateDigitVecCCW(digits, maxRes, INVALID_DIGIT);
        }
    }
}

} // namespace z7
} // namespace hexify
