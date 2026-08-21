# tests/testthat/test-namespace-hygiene.R
# Guards on the generated Rcpp bindings.
#
# Rcpp::compileAttributes() writes one R function per [[Rcpp::export]] into
# R/RcppExports.R. That file sorts before most of R/, so a binding named after
# a package function would be the definition R keeps or drops depending on
# collation order alone. Prefixing every entry point cpp_ keeps the two name
# spaces apart.

rcpp_binding_names <- function() {
  dll <- .getNamespaceInfo(asNamespace("hexify"), "DLLs")[["hexify"]]
  sub("^_hexify_", "", names(getDLLRegisteredRoutines(dll)$.Call))
}

test_that("every compiled entry point is cpp_ prefixed", {
  bindings <- rcpp_binding_names()
  expect_gt(length(bindings), 0)
  expect_equal(bindings[!startsWith(bindings, "cpp_")], character(0))
})

test_that("no compiled entry point shadows a package function", {
  clashes <- intersect(rcpp_binding_names(), getNamespaceExports("hexify"))
  expect_equal(clashes, character(0))
})

test_that("no compiled entry point shadows an internal R function", {
  ns <- asNamespace("hexify")
  r_sources <- setdiff(ls(ns, all.names = TRUE), rcpp_binding_names())
  expect_equal(intersect(rcpp_binding_names(), r_sources), character(0))
})
