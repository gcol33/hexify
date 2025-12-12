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
#include "cell_index.h"
#include "index_z3.h"
#include "index_z7.h"

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
  long long i, j;
  hexify::hex_quantize_ap3(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(fwd.face, i, j, resolution, 3, idx_type);
}

// [[Rcpp::export]]
std::string cpp_lonlat_to_index_ap4(double lon_deg, double lat_deg, int resolution,
                                     std::string index_type = "auto") {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i, j;
  hexify::hex_quantize_ap4(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(fwd.face, i, j, resolution, 4, idx_type);
}

// [[Rcpp::export]]
std::string cpp_lonlat_to_index_ap7(double lon_deg, double lat_deg, int resolution,
                                     std::string index_type = "auto") {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i, j;
  hexify::hex_quantize_ap7(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  hexify::IndexType idx_type = parse_index_type(index_type);
  return hexify::cell_to_index(fwd.face, i, j, resolution, 7, idx_type);
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

  double cx, cy;
  if (aperture == 3) {
    hexify::hex_center_ap3(i, j, resolution, cx, cy);
  } else if (aperture == 4) {
    hexify::hex_center_ap4(i, j, resolution, cx, cy);
  } else if (aperture == 7) {
    hexify::hex_center_ap7(i, j, resolution, cx, cy);
  } else {
    Rcpp::stop("Invalid aperture");
  }

  auto ll = hexify::face_xy_to_ll(cx, cy, face);

  return Rcpp::NumericVector::create(
    _["lon"] = ll.first,
    _["lat"] = ll.second
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

// ============================================================================
// Z3 Helpers
// ============================================================================

// [[Rcpp::export]]
List cpp_hex_index_z3_quantize_digits(double tx, double ty, int eff_res,
                                      double center_thr = 0.4,
                                      LogicalVector flip_classes = LogicalVector::create()) {
  // Placeholder implementation for Z3 quantization
  std::vector<int> digits;

  for (int i = 0; i < eff_res; i++) {
    digits.push_back(0);
  }

  return List::create(
    Named("digits") = digits,
    Named("tx") = tx,
    Named("ty") = ty
  );
}

// [[Rcpp::export]]
List cpp_hex_index_z3_center(IntegerVector d,
                             LogicalVector flip_classes = LogicalVector::create()) {
  double cx = 0.0;
  double cy = 0.0;

  return List::create(
    Named("cx") = cx,
    Named("cy") = cy
  );
}

// [[Rcpp::export]]
List cpp_hex_index_z3_corners(IntegerVector digs,
                              LogicalVector flip_classes = LogicalVector::create(),
                              double hex_radius = 1.0) {
  std::vector<double> x_coords;
  std::vector<double> y_coords;

  for (int i = 0; i < 6; i++) {
    double angle = hexify::kPiOver3 * i;
    x_coords.push_back(hex_radius * std::cos(angle));
    y_coords.push_back(hex_radius * std::sin(angle));
  }

  return List::create(
    Named("x") = x_coords,
    Named("y") = y_coords
  );
}
