// tri_to_q2di.h - Convert PROJTRI (triangle coords) to Q2DI (quad integer coords)
//
// Coordinate transformations for ISEA DGGS grids.
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_TRI_TO_Q2DI_H
#define HEXIFY_TRI_TO_Q2DI_H

namespace hexify {

// Convert from PROJTRI (tnum, tx, ty) to Q2DD (quad, qx, qy)
// This applies triTable rotation and translation
void projtri_to_q2dd(int tnum, double tx, double ty,
                     int& out_quad, double& out_qx, double& out_qy);

// Convert from Q2DD (quad, qx, qy) to Q2DI (quad, i, j)
// This quantizes continuous coords to integer cell indices
// Note: aperture and resolution affect the scaling
void q2dd_to_q2di(int quad, double qx, double qy,
                  int aperture, int resolution,
                  int& out_quad, long long& out_i, long long& out_j);

// Full pipeline: PROJTRI → Q2DI
void projtri_to_q2di(int tnum, double tx, double ty,
                     int aperture, int resolution,
                     int& out_quad, long long& out_i, long long& out_j);

// Inverse: Q2DI → Q2DD (for computing cell centers)
void q2di_to_q2dd(int quad, long long i, long long j,
                  int aperture, int resolution,
                  double& out_qx, double& out_qy);

// Inverse: Q2DD → PROJTRI (throws on invalid region)
void q2dd_to_projtri(int quad, double qx, double qy,
                     int& out_tnum, double& out_tx, double& out_ty);

// Inverse: Q2DD → PROJTRI (returns false on invalid region)
bool try_q2dd_to_projtri(int quad, double qx, double qy,
                         int& out_tnum, double& out_tx, double& out_ty);

// Get maxI/maxJ for a given aperture and resolution
long long get_max_ij(int aperture, int resolution);

// Edge handling: map edge cells to adjacent quads
// Returns true if coord was on edge and was adjusted
bool handle_edge_overflow(int& quadNum, long long& i, long long& j,
                          int aperture, int resolution);

} // namespace hexify

#endif // HEXIFY_TRI_TO_Q2DI_H
