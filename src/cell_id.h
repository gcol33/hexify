#pragma once
// cell_id.h - ISEA cell IDs at the R boundary
//
// An ISEA cell ID is a signed 64-bit integer, 1 to the grid's cell count. It
// crosses into R as bit64's integer64: the 64 bits of the integer stored in
// the slot of a double, with class "integer64", NA the smallest int64. A
// double holds whole numbers exactly only below 2^53, so IDs are never
// converted to one.

#include <Rcpp.h>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace hexify {

constexpr int64_t kCellIdNA = std::numeric_limits<int64_t>::min();
constexpr int64_t kCellIdMax = std::numeric_limits<int64_t>::max();

// The ID stored in an integer64 slot
inline int64_t cell_id_get(double slot) {
  int64_t v;
  std::memcpy(&v, &slot, sizeof v);
  return v;
}

// The integer64 slot holding an ID
inline double cell_id_slot(int64_t v) {
  double d;
  std::memcpy(&d, &v, sizeof d);
  return d;
}

// Stops unless `x` is an integer64 vector. A plain double reaching here would
// be read as the bits of an integer, so every entry point taking IDs checks.
inline void require_cell_ids(const Rcpp::NumericVector& x) {
  if (!Rf_inherits(x, "integer64")) {
    Rcpp::stop("hexify internal error: cell IDs must reach C++ as integer64");
  }
}

// An integer64 vector of n IDs, all NA
inline Rcpp::NumericVector cell_id_na(R_xlen_t n) {
  Rcpp::NumericVector out(n, cell_id_slot(kCellIdNA));
  out.attr("class") = "integer64";
  return out;
}

// An integer64 vector holding `ids`
inline Rcpp::NumericVector cell_id_vector(const std::vector<int64_t>& ids) {
  Rcpp::NumericVector out(ids.size());
  for (size_t k = 0; k < ids.size(); k++) out[k] = cell_id_slot(ids[k]);
  out.attr("class") = "integer64";
  return out;
}

} // namespace hexify
