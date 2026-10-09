// rcpp_hex9.cpp - Hex9 hierarchy and labels for R
//
// A Hex9 cell ID is its address (hex9.h), and the hierarchy and labels are
// read from the lattice alone: none of these entry points depends on where
// the octahedron sits or how its faces are projected.

#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "cell_id.h"
#include "hex9.h"
#include "hex9_warp.h"

using namespace Rcpp;

namespace {

void check_level(int resolution) {
  if (resolution < 0 || resolution > hexify::hex9::kMaxLevel) {
    stop("Hex9 resolution must be between 0 and %d", hexify::hex9::kMaxLevel);
  }
}

// The centre of the cell an ID names; stops on an ID naming none
hexify::hex9::OctPoint cell_centre(int64_t id, int resolution) {
  hexify::hex9::OctPoint c;
  if (!hexify::hex9::decode(id, resolution, c)) {
    stop("cell_id must be a whole number in [1, %lld] for Hex9 resolution %d",
         static_cast<long long>(hexify::hex9::level_cells(resolution)), resolution);
  }
  return c;
}

const char kDigit[] = "0123456789ab";

} // anon

// The canonical parent of each cell, one resolution up
// [[Rcpp::export]]
NumericVector cpp_hex9_parent(NumericVector cell_id, int resolution) {
  check_level(resolution);
  if (resolution < 1) stop("Hex9 cells of resolution 0 have no parent");
  hexify::require_cell_ids(cell_id);
  NumericVector out = hexify::cell_id_na(cell_id.size());
  for (R_xlen_t k = 0; k < cell_id.size(); k++) {
    const int64_t id = hexify::cell_id_get(cell_id[k]);
    if (id == hexify::kCellIdNA) continue;
    const hexify::hex9::OctPoint p = hexify::hex9::parent(cell_centre(id, resolution), resolution);
    out[k] = hexify::cell_id_slot(hexify::hex9::encode(p, resolution - 1));
  }
  return out;
}

// The nine canonical children of each cell, one resolution down, sorted
// [[Rcpp::export]]
List cpp_hex9_children(NumericVector cell_id, int resolution) {
  check_level(resolution + 1);
  hexify::require_cell_ids(cell_id);
  List out(cell_id.size());
  hexify::hex9::OctPoint kids[9];
  for (R_xlen_t k = 0; k < cell_id.size(); k++) {
    const int64_t id = hexify::cell_id_get(cell_id[k]);
    if (id == hexify::kCellIdNA) {
      out[k] = hexify::cell_id_na(0);
      continue;
    }
    hexify::hex9::children(cell_centre(id, resolution), resolution, kids);
    std::vector<int64_t> ids(9);
    for (int c = 0; c < 9; c++) ids[c] = hexify::hex9::encode(kids[c], resolution + 1);
    std::sort(ids.begin(), ids.end());
    out[k] = hexify::cell_id_vector(ids);
  }
  return out;
}

// Each cell's libhex9 label, its digits and key tail, "<digits>.<tail>"
// [[Rcpp::export]]
CharacterVector cpp_hex9_label(NumericVector cell_id, int resolution) {
  check_level(resolution);
  hexify::require_cell_ids(cell_id);
  CharacterVector out(cell_id.size());
  int digits[hexify::hex9::kMaxLevel + 1];
  for (R_xlen_t k = 0; k < cell_id.size(); k++) {
    const int64_t id = hexify::cell_id_get(cell_id[k]);
    if (id == hexify::kCellIdNA) {
      out[k] = NA_STRING;
      continue;
    }
    const hexify::hex9::OctPoint c = cell_centre(id, resolution);
    hexify::hex9::id_digits(id, resolution, digits);
    std::string s;
    for (int l = 0; l <= resolution; l++) s += kDigit[digits[l]];
    s += '.';
    s += static_cast<char>('0' + hexify::hex9::key_tail(c, resolution));
    out[k] = s;
  }
  return out;
}

// Labels back to cells: the resolution each names, its length less one, and
// the cell ID. A label may leave out its tail; one that names no cell, or
// whose tail is not the cell's, comes back NA.
// [[Rcpp::export]]
List cpp_hex9_parse_label(CharacterVector label) {
  const R_xlen_t n = label.size();
  IntegerVector res(n, NA_INTEGER);
  NumericVector ids = hexify::cell_id_na(n);
  for (R_xlen_t k = 0; k < n; k++) {
    if (CharacterVector::is_na(label[k])) continue;
    const std::string s = Rcpp::as<std::string>(label[k]);
    const size_t dot = s.find('.');
    const std::string body = s.substr(0, dot);
    const int level = static_cast<int>(body.size()) - 1;
    if (level < 0 || level > hexify::hex9::kMaxLevel) continue;
    int64_t id = 0;
    bool ok = true;
    for (int l = 0; l <= level && ok; l++) {
      const char* at = std::find(kDigit, kDigit + 12, body[l]);
      const int d = static_cast<int>(at - kDigit);
      if (d >= 12 || (l > 0 && d >= 9)) ok = false;
      id = id * 9 + d;
    }
    if (!ok) continue;
    id += 1;
    hexify::hex9::OctPoint c;
    if (!hexify::hex9::decode(id, level, c)) continue;
    if (dot != std::string::npos) {
      const std::string tail = s.substr(dot + 1);
      if (tail.size() != 1 || tail[0] - '0' != hexify::hex9::key_tail(c, level)) continue;
    }
    res[k] = level;
    ids[k] = hexify::cell_id_slot(id);
  }
  return List::create(_["resolution"] = res, _["cell_id"] = ids);
}

// The cell of each point (x, y, z) of the octahedron |x| + |y| + |z| = 1 in
// the solid's standard frame, with no projection: the combinatorial grid
// alone, as libhex9 reads its own face charts.
// [[Rcpp::export]]
NumericVector cpp_hex9_octahedron_cell(NumericVector x, NumericVector y,
                                       NumericVector z, int resolution) {
  check_level(resolution);
  const R_xlen_t n = x.size();
  if (y.size() != n || z.size() != n) stop("x, y and z must be the same length");
  NumericVector out = hexify::cell_id_na(n);
  for (R_xlen_t k = 0; k < n; k++) {
    const double p[3] = {x[k], y[k], z[k]};
    const double s1 = std::fabs(p[0]) + std::fabs(p[1]) + std::fabs(p[2]);
    if (!(s1 > 0.0) || !std::isfinite(s1)) continue;
    double w[3];
    int s[3];
    for (int a = 0; a < 3; a++) {
      w[a] = std::fabs(p[a]) / s1;
      s[a] = p[a] < 0.0 ? -1 : 1;
    }
    const hexify::hex9::OctPoint c = hexify::hex9::locate(w, s, resolution);
    out[k] = hexify::cell_id_slot(hexify::hex9::encode(c, resolution));
  }
  return out;
}

// Load Hex9's warp field from the bytes of libhex9's .h9warp file
// [[Rcpp::export]]
void cpp_hex9_warp_load(RawVector bytes) {
  std::string err;
  if (!hexify::hex9::warp_load(RAW(bytes), static_cast<std::size_t>(bytes.size()), err)) {
    stop(err);
  }
}

// Whether Hex9's warp field is loaded
// [[Rcpp::export]]
bool cpp_hex9_warp_ready() {
  return hexify::hex9::warp_ready();
}

// Whether the Hex9 digit tables give every cell of a level its own digit
// string at every level (hexify::hex9::digit_strings_unique())
// [[Rcpp::export]]
bool cpp_hex9_digit_strings_unique() {
  return hexify::hex9::digit_strings_unique();
}
