// rcpp_index.cpp
// Rcpp bindings for cell indexing and hierarchy functions
//
// This file provides the R interface for:
// - Cell to index conversion
// - Index to cell conversion
// - Parent/child hierarchy navigation
// - Index comparison and validation
// - Z7 canonical form

#include <Rcpp.h>
#include "constants.h"
#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "aperture.h"
#include "coordinate_transforms.h"
#include "cell_index.h"
#include "index_z3.h"
#include "index_z7.h"
#include <cstdlib>

using namespace Rcpp;

// Helper to convert R string to IndexType enum
static hexify::IndexType parse_index_type(const std::string& type_str) {
  if (type_str == "auto" || type_str == "AUTO") return hexify::IndexType::AUTO;
  if (type_str == "zorder" || type_str == "ZORDER") return hexify::IndexType::ZORDER;
  if (type_str == "z3" || type_str == "Z3") return hexify::IndexType::Z3;
  if (type_str == "z7" || type_str == "Z7") return hexify::IndexType::Z7;
  Rcpp::stop("Invalid index_type. Must be 'auto', 'zorder', 'z3', or 'z7'");
}

// ============================================================================
// Cell to Index Conversion
// ============================================================================

// [[Rcpp::export]]
std::string cpp_cell_to_index(int face, double i, double j,
                               int resolution, int aperture,
                               std::string index_type = "auto") {
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(face,
                                static_cast<long long>(i),
                                static_cast<long long>(j),
                                resolution, aperture, idx_type);
}

// [[Rcpp::export]]
Rcpp::List cpp_index_to_cell(std::string index, int aperture,
                              std::string index_type = "auto") {
  int face, resolution;
  long long i, j;

  hexify::IndexType idx_type = parse_index_type(index_type);
  hexify::index_to_cell(index, aperture, idx_type, face, i, j, resolution);

  return Rcpp::List::create(
    Rcpp::Named("face") = face,
    Rcpp::Named("i") = static_cast<double>(i),
    Rcpp::Named("j") = static_cast<double>(j),
    Rcpp::Named("resolution") = resolution
  );
}

// ============================================================================
// Hierarchy Navigation
// ============================================================================

// [[Rcpp::export]]
std::string cpp_get_parent_index(std::string index, int aperture,
                                  std::string index_type = "auto") {
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::get_parent_index(index, aperture, idx_type);
}

// [[Rcpp::export]]
Rcpp::StringVector cpp_get_children_indices(std::string index, int aperture,
                                             std::string index_type = "auto") {
  hexify::IndexType idx_type = parse_index_type(index_type);
  std::vector<std::string> children = hexify::get_children_indices(index, aperture, idx_type);
  return Rcpp::wrap(children);
}

// [[Rcpp::export]]
int cpp_get_index_resolution(std::string index, int aperture,
                              std::string index_type = "auto") {
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::get_index_resolution(index, aperture, idx_type);
}

// ============================================================================
// Index Comparison and Validation
// ============================================================================

// [[Rcpp::export]]
int cpp_compare_indices(std::string idx1, std::string idx2) {
  return hexify::compare_indices(idx1, idx2);
}

// [[Rcpp::export]]
bool cpp_is_valid_index_type(int aperture, std::string index_type) {
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::is_valid_index_type(aperture, idx_type);
}

// [[Rcpp::export]]
std::string cpp_get_default_index_type(int aperture) {
  hexify::IndexType idx_type = hexify::get_default_index_type(aperture);
  switch(idx_type) {
    case hexify::IndexType::ZORDER: return "zorder";
    case hexify::IndexType::Z3: return "z3";
    case hexify::IndexType::Z7: return "z7";
    default: return "auto";
  }
}

// ============================================================================
// Lon/Lat to Index Conversion (Aperture-specific)
// ============================================================================

// [[Rcpp::export]]
std::string cpp_lonlat_to_index_ap3(double lon_deg, double lat_deg, int resolution,
                                     std::string index_type = "auto") {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  // Quantizing directly against the raw triangle-frame coordinates (as this
  // used to do) yields (i,j) centered on the triangle's local origin, which
  // can be negative -- a different, ungrounded coordinate frame from what
  // cell_to_index()/z3::encode() expect (non-negative offsets from the
  // face's corner, as produced by cpp_lonlat_to_cell() via the same fold
  // used for aperture 4/7 below). Route through icosa_tri_to_quad_ij() first,
  // matching cpp_lonlat_to_index_ap4/ap7() and cpp_lonlat_to_cell().
  int face;
  long long i, j;
  hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                3, resolution, face, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(face, i, j, resolution, 3, idx_type);
}

// [[Rcpp::export]]
std::string cpp_lonlat_to_index_ap4(double lon_deg, double lat_deg, int resolution,
                                     std::string index_type = "auto") {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  // cell_to_index() only accepts face 0-11 for aperture 4/7 (a "quad", not the
  // raw 0-19 icosahedron triangle), so the triangle must first be folded into
  // its quad via icosa_tri_to_quad_ij() -- the same step cpp_lonlat_to_cell()
  // performs. Quantizing directly against the triangle-frame coordinates and
  // passing fwd.face through unfolded (as this used to do) throws for any
  // point on a triangle numbered 12-19.
  int quad;
  long long i, j;
  hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                4, resolution, quad, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(quad, i, j, resolution, 4, idx_type);
}

// [[Rcpp::export]]
std::string cpp_lonlat_to_index_ap7(double lon_deg, double lat_deg, int resolution,
                                     std::string index_type = "auto") {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  int quad;
  long long i, j;
  hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                7, resolution, quad, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(quad, i, j, resolution, 7, idx_type);
}

// [[Rcpp::export]]
std::string cpp_lonlat_to_index(double lon_deg, double lat_deg,
                                 int resolution, int aperture,
                                 std::string index_type = "auto") {
  if (aperture == 3) {
    return cpp_lonlat_to_index_ap3(lon_deg, lat_deg, resolution, index_type);
  } else if (aperture == 4) {
    return cpp_lonlat_to_index_ap4(lon_deg, lat_deg, resolution, index_type);
  } else if (aperture == 7) {
    return cpp_lonlat_to_index_ap7(lon_deg, lat_deg, resolution, index_type);
  } else {
    Rcpp::stop("Invalid aperture. Must be 3, 4, or 7");
  }
}

// ============================================================================
// Index to Lon/Lat Conversion
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector cpp_index_to_lonlat(std::string index, int aperture,
                                         std::string index_type = "auto") {
  int face, resolution;
  long long i, j;

  hexify::IndexType idx_type = parse_index_type(index_type);
  hexify::index_to_cell(index, aperture, idx_type, face, i, j, resolution);

  // `face` decoded from the index is a quad (0-11 for aperture 3 too, now
  // that cpp_lonlat_to_index_ap3() folds through icosa_tri_to_quad_ij()
  // like ap4/ap7 do), not a raw icosahedron triangle face -- it must be
  // folded back to the triangle frame via
  // quad_ij_to_xy()/quad_xy_to_icosa_tri() before face_xy_to_ll(),
  // mirroring the forward fold done in cpp_lonlat_to_index_ap3/ap4/ap7().
  double quad_x, quad_y;
  hexify::quad_ij_to_xy(face, i, j, aperture, resolution, quad_x, quad_y);

  int tri_face;
  double tri_x, tri_y;
  hexify::quad_xy_to_icosa_tri(face, quad_x, quad_y, tri_face, tri_x, tri_y);

  auto ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
  double lon = ll.first;
  double lat = ll.second;

  return Rcpp::NumericVector::create(
    _["lon"] = lon,
    _["lat"] = lat
  );
}

// ============================================================================
// Legacy exports (for backwards compatibility)
// ============================================================================

// [[Rcpp::export]]
std::string cell_to_index(int face, Rcpp::NumericVector i, Rcpp::NumericVector j,
                          int resolution, int aperture, std::string index_type = "AUTO") {
  hexify::IndexType type = parse_index_type(index_type);

  long long i_val = static_cast<long long>(i[0]);
  long long j_val = static_cast<long long>(j[0]);

  return hexify::cell_to_index(face, i_val, j_val, resolution, aperture, type);
}

// [[Rcpp::export]]
Rcpp::List index_to_cell(std::string index, int aperture, std::string index_type = "AUTO") {
  hexify::IndexType type = parse_index_type(index_type);

  int face, resolution;
  long long i, j;

  hexify::index_to_cell(index, aperture, type, face, i, j, resolution);

  return Rcpp::List::create(
    Rcpp::Named("face") = face,
    Rcpp::Named("i") = static_cast<double>(i),
    Rcpp::Named("j") = static_cast<double>(j),
    Rcpp::Named("resolution") = resolution
  );
}

// [[Rcpp::export]]
std::string get_parent_index(std::string index, int aperture, std::string index_type = "AUTO") {
  hexify::IndexType type = parse_index_type(index_type);
  return hexify::get_parent_index(index, aperture, type);
}

// [[Rcpp::export]]
std::vector<std::string> get_children_indices(std::string index, int aperture,
                                               std::string index_type = "AUTO") {
  hexify::IndexType type = parse_index_type(index_type);
  return hexify::get_children_indices(index, aperture, type);
}

// [[Rcpp::export]]
int get_index_resolution(std::string index, int aperture, std::string index_type = "AUTO") {
  hexify::IndexType type = parse_index_type(index_type);
  return hexify::get_index_resolution(index, aperture, type);
}

// ============================================================================
// Z7 Canonical Form
// ============================================================================

// [[Rcpp::export]]
std::string cpp_z7_canonical_form(std::string index, int max_iterations = 128) {
  try {
    return hexify::z7::canonical_form(index, max_iterations);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in z7_canonical_form: %s", e.what());
  }
}

