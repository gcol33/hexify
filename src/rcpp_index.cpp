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

  double lon, lat;
  if (aperture == 3) {
    double cx, cy;
    hexify::hex_center_ap3(i, j, resolution, cx, cy);
    auto ll = hexify::face_xy_to_ll(cx, cy, face);
    lon = ll.first;
    lat = ll.second;
  } else {
    // For aperture 4/7, `face` decoded from the index is a quad (0-11), not
    // a raw icosahedron triangle face -- it must be folded back to the
    // triangle frame via quad_ij_to_xy()/quad_xy_to_icosa_tri() before
    // face_xy_to_ll(), mirroring the forward fold done in
    // cpp_lonlat_to_index_ap4/ap7().
    double quad_x, quad_y;
    hexify::quad_ij_to_xy(face, i, j, aperture, resolution, quad_x, quad_y);

    int tri_face;
    double tri_x, tri_y;
    hexify::quad_xy_to_icosa_tri(face, quad_x, quad_y, tri_face, tri_x, tri_y);

    auto ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
    lon = ll.first;
    lat = ll.second;
  }

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

// ============================================================================
// Z3 Helpers
//
// These bind hexify_assign()'s per-point Z3 workflow (quantize -> digits ->
// center/corners) onto the same quantization primitives and Z3 digit table
// used everywhere else Z3 cell IDs are produced (hexify::hex_quantize_ap3/
// hex_center_ap3/hex_corners_ap3 and hexify::z3::encode/decode), rather than
// re-deriving the geometry. z3::encode(i, j, resolution) always returns a
// string of exactly `resolution` characters ('0'-'2'), so it maps 1:1 onto a
// "digits" vector of the same length -- no separate digit scheme is needed.
//
// `flip_classes`/`center_thr` are accepted for API compatibility with the
// documented `match_dggrid_parity` argument but are not yet wired to a
// concrete effect: no DGGRID ISEA3H parity convention distinct from the
// Class I/II alternation already implemented in hex_quantize_ap3/z3::encode
// has been specified. See the tracking issue for this gap.
// ============================================================================

namespace {

std::string digits_to_z3_string(const IntegerVector& d) {
  std::string s;
  s.reserve(d.size());
  for (int idx = 0; idx < d.size(); ++idx) {
    int v = d[idx];
    if (v < 0 || v > 2) {
      Rcpp::stop("hex_index_z3: digit at position %d must be 0, 1, or 2 (got %d)", idx + 1, v);
    }
    s += static_cast<char>('0' + v);
  }
  return s;
}

} // anonymous namespace

// [[Rcpp::export]]
List cpp_hex_index_z3_quantize_digits(double tx, double ty, int eff_res,
                                      double center_thr = 0.4,
                                      LogicalVector flip_classes = LogicalVector::create()) {
  long long i, j;
  hexify::hex_quantize_ap3(tx, ty, eff_res, i, j);

  std::string z3_str = hexify::z3::encode(i, j, eff_res);

  std::vector<int> digits;
  digits.reserve(z3_str.size());
  for (char c : z3_str) {
    digits.push_back(c - '0');
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
  std::string z3_str = digits_to_z3_string(d);
  int resolution = static_cast<int>(d.size());

  long long i, j;
  hexify::z3::decode(z3_str, resolution, i, j);

  double cx, cy;
  hexify::hex_center_ap3(i, j, resolution, cx, cy);

  return List::create(
    Named("cx") = cx,
    Named("cy") = cy
  );
}

// [[Rcpp::export]]
List cpp_hex_index_z3_corners(IntegerVector digs,
                              LogicalVector flip_classes = LogicalVector::create(),
                              double hex_radius = 1.0) {
  std::string z3_str = digits_to_z3_string(digs);
  int resolution = static_cast<int>(digs.size());

  long long i, j;
  hexify::z3::decode(z3_str, resolution, i, j);

  double out_x[6], out_y[6];
  hexify::hex_corners_ap3(i, j, resolution, hex_radius, out_x, out_y);

  std::vector<double> x_coords(out_x, out_x + 6);
  std::vector<double> y_coords(out_y, out_y + 6);

  return List::create(
    Named("x") = x_coords,
    Named("y") = y_coords
  );
}
