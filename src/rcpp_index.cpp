// rcpp_index.cpp
// Rcpp bindings for cell indexing and hierarchy functions
//
// This file provides the R interface for:
// - Cell to index conversion
// - Index to cell conversion
// - Parent/child hierarchy navigation
// - Index comparison and validation
// - Z7 canonical form
//
// Each entry point takes vectors and loops here, so a caller holding a vector
// of cells or indices pays one call rather than one per element. A missing
// input carries through as a missing result rather than reaching the encoder.

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

// Two vector arguments of one call either run in step, or one of them is a
// single value read for every element.
static R_xlen_t paired_length(R_xlen_t n, R_xlen_t m, const char* fn,
                              const char* a, const char* b) {
  if (n == m) return n;
  if (n == 1) return m;
  if (m == 1) return n;
  Rcpp::stop("%s: `%s` (%d) and `%s` (%d) must be the same length, or one of "
             "them length 1", fn, a, static_cast<int>(n), b, static_cast<int>(m));
}

static void require_same_length(R_xlen_t n, R_xlen_t m, const char* fn,
                                const char* a, const char* b) {
  if (n != m) {
    Rcpp::stop("%s: `%s` (%d) and `%s` (%d) must be the same length",
               fn, a, static_cast<int>(n), b, static_cast<int>(m));
  }
}

// ============================================================================
// Cell to Index Conversion
// ============================================================================

// [[Rcpp::export]]
CharacterVector cpp_cell_to_index(IntegerVector face, NumericVector i,
                                   NumericVector j, int resolution,
                                   int aperture,
                                   std::string index_type = "auto") {
  R_xlen_t n = face.size();
  require_same_length(n, i.size(), "cpp_cell_to_index", "face", "i");
  require_same_length(n, j.size(), "cpp_cell_to_index", "face", "j");

  hexify::IndexType idx_type = parse_index_type(index_type);
  CharacterVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (IntegerVector::is_na(face[k]) || NumericVector::is_na(i[k]) ||
        NumericVector::is_na(j[k])) {
      out[k] = NA_STRING;
      continue;
    }
    out[k] = hexify::cell_to_index(face[k],
                                    static_cast<long long>(i[k]),
                                    static_cast<long long>(j[k]),
                                    resolution, aperture, idx_type);
  }

  return out;
}

// [[Rcpp::export]]
DataFrame cpp_index_to_cell(CharacterVector index, int aperture,
                             std::string index_type = "auto") {
  R_xlen_t n = index.size();
  hexify::IndexType idx_type = parse_index_type(index_type);

  IntegerVector out_face(n);
  NumericVector out_i(n);
  NumericVector out_j(n);
  IntegerVector out_res(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out_face[k] = NA_INTEGER;
      out_i[k] = NA_REAL;
      out_j[k] = NA_REAL;
      out_res[k] = NA_INTEGER;
      continue;
    }

    int face, resolution;
    long long i, j;
    hexify::index_to_cell(Rcpp::as<std::string>(index[k]), aperture, idx_type,
                          face, i, j, resolution);

    out_face[k] = face;
    out_i[k] = static_cast<double>(i);
    out_j[k] = static_cast<double>(j);
    out_res[k] = resolution;
  }

  return DataFrame::create(
    _["face"] = out_face,
    _["i"] = out_i,
    _["j"] = out_j,
    _["resolution"] = out_res
  );
}

// ============================================================================
// Hierarchy Navigation
// ============================================================================

// [[Rcpp::export]]
CharacterVector cpp_get_parent_index(CharacterVector index, int aperture,
                                      std::string index_type = "auto") {
  R_xlen_t n = index.size();
  hexify::IndexType idx_type = parse_index_type(index_type);
  CharacterVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out[k] = NA_STRING;
      continue;
    }
    out[k] = hexify::get_parent_index(Rcpp::as<std::string>(index[k]),
                                       aperture, idx_type);
  }

  return out;
}

// [[Rcpp::export]]
List cpp_get_children_indices(CharacterVector index, int aperture,
                               std::string index_type = "auto") {
  R_xlen_t n = index.size();
  hexify::IndexType idx_type = parse_index_type(index_type);
  List out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out[k] = CharacterVector(0);
      continue;
    }
    std::vector<std::string> children =
      hexify::get_children_indices(Rcpp::as<std::string>(index[k]),
                                    aperture, idx_type);
    out[k] = Rcpp::wrap(children);
  }

  return out;
}

// [[Rcpp::export]]
IntegerVector cpp_get_index_resolution(CharacterVector index, int aperture,
                                        std::string index_type = "auto") {
  R_xlen_t n = index.size();
  hexify::IndexType idx_type = parse_index_type(index_type);
  IntegerVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out[k] = NA_INTEGER;
      continue;
    }
    out[k] = hexify::get_index_resolution(Rcpp::as<std::string>(index[k]),
                                           aperture, idx_type);
  }

  return out;
}

// ============================================================================
// Index Comparison and Validation
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_compare_indices(CharacterVector idx1, CharacterVector idx2) {
  R_xlen_t n1 = idx1.size(), n2 = idx2.size();
  R_xlen_t n = paired_length(n1, n2, "cpp_compare_indices", "idx1", "idx2");
  IntegerVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    R_xlen_t k1 = (n1 == 1) ? 0 : k;
    R_xlen_t k2 = (n2 == 1) ? 0 : k;
    if (CharacterVector::is_na(idx1[k1]) || CharacterVector::is_na(idx2[k2])) {
      out[k] = NA_INTEGER;
      continue;
    }
    out[k] = hexify::compare_indices(Rcpp::as<std::string>(idx1[k1]),
                                      Rcpp::as<std::string>(idx2[k2]));
  }

  return out;
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
// Lon/Lat to Index Conversion
// ============================================================================

// cell_to_index()/z3::encode() read (i,j) as non-negative offsets from the
// corner of a quad, and for aperture 4 and 7 they take a quad (0-11) rather
// than a raw icosahedron triangle (0-19). Quantizing the triangle-frame
// coordinates directly gives (i,j) centered on the triangle's local origin,
// which can be negative, and leaves a point on triangle 12-19 carrying a face
// number the encoder rejects. icosa_tri_to_quad_ij() performs the fold, the
// same step cpp_lonlat_to_cell() takes.
static std::string lonlat_to_index_one(double lon_deg, double lat_deg,
                                        int resolution, int aperture,
                                        hexify::IndexType idx_type) {
  hexify::ProjectionResult fwd = hexify::snyder_forward(lon_deg, lat_deg);
  int quad;
  long long i, j;
  hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x,
                                fwd.icosa_triangle_y,
                                aperture, resolution, quad, i, j);
  return hexify::cell_to_index(quad, i, j, resolution, aperture, idx_type);
}

static CharacterVector lonlat_to_index_vec(const NumericVector& lon,
                                            const NumericVector& lat,
                                            int resolution, int aperture,
                                            const std::string& index_type,
                                            const char* fn) {
  R_xlen_t n = lon.size();
  require_same_length(n, lat.size(), fn, "lon", "lat");

  hexify::IndexType idx_type = parse_index_type(index_type);
  CharacterVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (NumericVector::is_na(lon[k]) || NumericVector::is_na(lat[k])) {
      out[k] = NA_STRING;
      continue;
    }
    out[k] = lonlat_to_index_one(lon[k], lat[k], resolution, aperture, idx_type);
  }

  return out;
}

// [[Rcpp::export]]
CharacterVector cpp_lonlat_to_index_ap3(NumericVector lon_deg,
                                         NumericVector lat_deg, int resolution,
                                         std::string index_type = "auto") {
  return lonlat_to_index_vec(lon_deg, lat_deg, resolution, 3, index_type,
                             "cpp_lonlat_to_index_ap3");
}

// [[Rcpp::export]]
CharacterVector cpp_lonlat_to_index_ap4(NumericVector lon_deg,
                                         NumericVector lat_deg, int resolution,
                                         std::string index_type = "auto") {
  return lonlat_to_index_vec(lon_deg, lat_deg, resolution, 4, index_type,
                             "cpp_lonlat_to_index_ap4");
}

// [[Rcpp::export]]
CharacterVector cpp_lonlat_to_index_ap7(NumericVector lon_deg,
                                         NumericVector lat_deg, int resolution,
                                         std::string index_type = "auto") {
  return lonlat_to_index_vec(lon_deg, lat_deg, resolution, 7, index_type,
                             "cpp_lonlat_to_index_ap7");
}

// [[Rcpp::export]]
CharacterVector cpp_lonlat_to_index(NumericVector lon_deg,
                                     NumericVector lat_deg,
                                     int resolution, int aperture,
                                     std::string index_type = "auto") {
  if (aperture != 3 && aperture != 4 && aperture != 7) {
    Rcpp::stop("Invalid aperture. Must be 3, 4, or 7");
  }
  return lonlat_to_index_vec(lon_deg, lat_deg, resolution, aperture, index_type,
                             "cpp_lonlat_to_index");
}

// ============================================================================
// Index to Lon/Lat Conversion
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_index_to_lonlat(CharacterVector index, int aperture,
                               std::string index_type = "auto") {
  R_xlen_t n = index.size();
  hexify::IndexType idx_type = parse_index_type(index_type);

  NumericVector out_lon(n);
  NumericVector out_lat(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out_lon[k] = NA_REAL;
      out_lat[k] = NA_REAL;
      continue;
    }

    int face, resolution;
    long long i, j;
    hexify::index_to_cell(Rcpp::as<std::string>(index[k]), aperture, idx_type,
                          face, i, j, resolution);

    // `face` decoded from the index is a quad (0-11 for aperture 3 too, since
    // the forward direction folds through icosa_tri_to_quad_ij()), not a raw
    // icosahedron triangle face -- it is folded back to the triangle frame
    // before face_xy_to_ll(), mirroring that forward fold.
    double quad_x, quad_y;
    hexify::quad_ij_to_xy(face, i, j, aperture, resolution, quad_x, quad_y);

    int tri_face;
    double tri_x, tri_y;
    hexify::quad_xy_to_icosa_tri(face, quad_x, quad_y, tri_face, tri_x, tri_y);

    std::pair<double, double> ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
    out_lon[k] = ll.first;
    out_lat[k] = ll.second;
  }

  return DataFrame::create(
    _["lon"] = out_lon,
    _["lat"] = out_lat
  );
}

// ============================================================================
// Z7 Canonical Form
// ============================================================================

// [[Rcpp::export]]
CharacterVector cpp_z7_canonical_form(CharacterVector index,
                                       int max_iterations = 128) {
  R_xlen_t n = index.size();
  CharacterVector out(n);

  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(index[k])) {
      out[k] = NA_STRING;
      continue;
    }
    try {
      out[k] = hexify::z7::canonical_form(Rcpp::as<std::string>(index[k]),
                                           max_iterations);
    } catch (const std::exception& e) {
      Rcpp::stop("Error in z7_canonical_form: %s", e.what());
    }
  }

  return out;
}
