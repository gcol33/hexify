// ijk_coordinates.h
// Cube Coordinates for Hexagonal Grids
//
// Implements the (i,j,k) cube coordinate system for hexagonal grids
// where the constraint i + j + k = 0 holds. This representation
// simplifies many hex operations including rotation and distance.
//
// The cube coordinate system embeds a 2D hex grid in 3D space along
// the plane x + y + z = 0. Each hex has 6 neighbors at unit distance.
//
// Aperture-7 operations:
// - upAp7/upAp7r: Coarsen by factor sqrt(7) with ~19.1° rotation
// - downAp7/downAp7r: Refine by factor sqrt(7)
// - The 'r' variants use opposite rotation direction
//
// References:
// - Red Blob Games "Hexagonal Grids" (redblobgames.com/grids/hexagons)
// - H3 Coordinate Systems (h3geo.org/docs/core-library/coordsystems)
// - Sahr et al. (2003) "Geodesic Discrete Global Grid Systems"
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once

#include <cstdint>

namespace hexify {
namespace z7 {

class IVec3D {
public:
  long long i_, j_, k_;
  
  enum Direction {
    CENTER_DIGIT = 0,   // (0, 0, 0)
    K_AXES_DIGIT = 1,   // (0, 0, 1)
    J_AXES_DIGIT = 2,   // (0, 1, 0)
    I_AXES_DIGIT = 3,   // (0, 1, 1)
    IK_AXES_DIGIT = 4,  // (1, 0, 0)
    IJ_AXES_DIGIT = 5,  // (1, 0, 1)
    JK_AXES_DIGIT = 6,  // (1, 1, 0)
    NUM_DIGITS = 7,
    INVALID_DIGIT = 7
  };
  
  static const int PENTAGON_SKIPPED_DIGIT_TYPE1 = 2;
  static const int PENTAGON_SKIPPED_DIGIT_TYPE2 = 5;
  
  IVec3D() : i_(0), j_(0), k_(0) {}
  IVec3D(long long i, long long j) : i_(i), j_(j), k_(-(i + j)) {}
  IVec3D(long long i, long long j, long long k) : i_(i), j_(j), k_(k) {}
  
  long long i() const { return i_; }
  long long j() const { return j_; }
  long long k() const { return k_; }
  
  void setI(long long val) { i_ = val; k_ = -(i_ + j_); }
  void setJ(long long val) { j_ = val; k_ = -(i_ + j_); }
  void setK(long long val) { k_ = val; }
  
  void ijkPlusNormalize();
  Direction unitIjkPlusToDigit() const;
  void upAp7();
  void upAp7r();
  void downAp7();
  void downAp7r();
  void neighbor(Direction digit);
  void ijkRotate60cw();
  void ijkRotate60ccw();
  
  IVec3D diffVec(const IVec3D& other) const {
    return IVec3D(i_ - other.i_, j_ - other.j_, k_ - other.k_);
  }
  
  static Direction rotate60ccw(Direction digit);
  static Direction rotate60cw(Direction digit);
  static void rotateDigitVecCCW(Direction* digits, int maxRes, Direction skipDigit);
  
  IVec3D operator+(const IVec3D& other) const {
    return IVec3D(i_ + other.i_, j_ + other.j_, k_ + other.k_);
  }
  
  IVec3D& operator+=(const IVec3D& other) {
    i_ += other.i_;
    j_ += other.j_;
    k_ += other.k_;
    return *this;
  }
  
  IVec3D operator*(long long scalar) const {
    return IVec3D(i_ * scalar, j_ * scalar, k_ * scalar);
  }
  
  IVec3D& operator*=(long long scalar) {
    i_ *= scalar;
    j_ *= scalar;
    k_ *= scalar;
    return *this;
  }
  
  bool operator==(const IVec3D& other) const {
    return i_ == other.i_ && j_ == other.j_ && k_ == other.k_;
  }
};

class IVec2D {
public:
  long long i_, j_;
  
  IVec2D(long long i, long long j) : i_(i), j_(j) {}
  explicit IVec2D(const IVec3D& ijk) : i_(ijk.i() - ijk.k()), j_(ijk.j() - ijk.k()) {}
  
  long long i() const { return i_; }
  long long j() const { return j_; }
  void setI(long long val) { i_ = val; }
  void setJ(long long val) { j_ = val; }
};

} // namespace z7
} // namespace hexify
