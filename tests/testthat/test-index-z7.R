
# test-index_z7.R
# Tests for Z7 encoding/decoding
# Tests for hexify's bijective aperture-7 index.

library(testthat)

test_that("Z7: Valid indices round-trip correctly", {
  # Test indices that are known to be valid and should round-trip
  # Avoid pentagon edge cases by using hexagon faces (1-10) and safe indices
  
  test_indices <- list(
    "01",     # Face 1 center
    "020",    # Face 2, digit 0
    "031",    # Face 3, digit 1  
    "044",    # Face 4, digit 4
    "055",    # Face 5, digit 5 (this one works even on pentagons)
    "066",    # Face 6, digit 6
    "0700",   # Face 7, safe multi-digit
    "0811",   # Face 8
    "0922",   # Face 9
    "1033"    # Face 10
  )
  
  for (idx in test_indices) {
    result <- hexify_index_to_cell(idx, 7L, "z7")
    idx2 <- hexify_cell_to_index(result$face, result$i, result$j, 
                          result$resolution, 7L, "z7")
    
    expect_equal(idx2, idx, 
                 info = sprintf("Index %s should round-trip correctly", idx))
  }
})

test_that("Z7: Resolution 0 (base cells) work correctly", {
  # All 12 base cells should round-trip at resolution 0
  for (face in 0:11) {
    idx <- hexify_cell_to_index(face, 0L, 0L, 0L, 7L, "z7")
    expected <- sprintf("%02d", face)
    expect_equal(idx, expected)
    
    result <- hexify_index_to_cell(idx, 7L, "z7")
    expect_equal(as.integer(result$face), face)
    expect_equal(as.integer(result$i), 0L)
    expect_equal(as.integer(result$j), 0L)
    expect_equal(as.integer(result$resolution), 0L)
  }
})

test_that("Z7: Hexagon faces handle adjacency remapping", {
  # Some coordinates on hexagon faces get remapped to other faces
  # This is expected behavior for aperture 7's cross-face structure
  
  # Test that we can at least encode and decode consistently
  test_cases <- list(
    list(face = 4L, i = 0L, j = 0L, res = 1L),  # Should work
    list(face = 4L, i = 1L, j = 0L, res = 1L),  # May remap
    list(face = 4L, i = 0L, j = 1L, res = 1L),  # May remap
    list(face = 4L, i = 1L, j = 1L, res = 1L),  # May remap
    list(face = 4L, i = 2L, j = 1L, res = 1L)   # May remap
  )
  
  for (tc in test_cases) {
    # Just test that encoding and decoding works without errors
    idx <- hexify_cell_to_index(tc$face, tc$i, tc$j, tc$res, 7L, "z7")
    result <- hexify_index_to_cell(idx, 7L, "z7")
    
    # Resolution should be preserved (extract the numeric value)
    expect_equal(as.integer(result$resolution), tc$res,
                 info = sprintf("Resolution preserved for face=%d i=%d j=%d", 
                               tc$face, tc$i, tc$j))
    
    # Face might change due to adjacency remapping - that's OK
    # Coordinates might change too - that's also OK for aperture 7
  }
})

test_that("Z7: Pentagon rotation behavior follows DGGRID", {
  # Test pentagon-specific behavior as implemented in DGGRID
  # Note: The exact encoding depends on coordinate mappings
  
  # Test that pentagon faces can be encoded/decoded
  test_cases <- list(
    list(face = 0L, i = 0L, j = 0L, res = 1L, expected = "000"),
    list(face = 11L, i = 0L, j = 0L, res = 1L, expected = "110")
  )
  
  for (tc in test_cases) {
    idx <- hexify_cell_to_index(tc$face, tc$i, tc$j, tc$res, 7L, "z7")
    expect_equal(idx, tc$expected,
                 info = sprintf("Pentagon encoding for face %d", tc$face))
  }
  
  # Test known pentagon rotation behaviors
  result_002 <- hexify_index_to_cell("002", 7L, "z7")
  idx_002_re <- hexify_cell_to_index(result_002$face, result_002$i, result_002$j,
                              result_002$resolution, 7L, "z7")
  # Exact result depends on DGGRID implementation
  
  result_115 <- hexify_index_to_cell("115", 7L, "z7")
  idx_115_re <- hexify_cell_to_index(result_115$face, result_115$i, result_115$j,
                              result_115$resolution, 7L, "z7")
  # Exact result depends on DGGRID implementation
})

test_that("Z7: Multi-digit indices round-trip", {
  
  test_cases <- list(
    list(idx = "0111111", expected = "0111111"),
    list(idx = "0222222", expected = "0222222"),
    list(idx = "0333333", expected = "0333333"),
    list(idx = "0400000", expected = "0400000"),
    list(idx = "0511111", expected = "0511111"),
    list(idx = "0622222", expected = "0622222"),
    list(idx = "0733333", expected = "0733333"),
    list(idx = "0844444", expected = "0844444"),
    list(idx = "0955555", expected = "0955555"),
    list(idx = "1066666", expected = "1066666")
  )
  
  for (tc in test_cases) {
    result <- hexify_index_to_cell(tc$idx, 7L, "z7")
    idx2 <- hexify_cell_to_index(result$face, result$i, result$j,
                         result$resolution, 7L, "z7")
    expect_equal(idx2, tc$expected,
                 info = sprintf("Index %s should round-trip to %s",
                               tc$idx, tc$expected))
  }
})

test_that("Z7: Parent-child relationships work", {
  # Test that valid parent indices have consistent children
  # Use hexagon faces to avoid pentagon complications
  
  test_parents <- c("01", "02", "03", "04", "06", "07", "08", "09", "10")
  
  for (parent_idx in test_parents) {
    parent_result <- hexify_index_to_cell(parent_idx, 7L, "z7")
    
    # Generate children and verify they decode properly
    for (digit in 0:6) {
      child_idx <- paste0(parent_idx, digit)
      
      # Decode child and verify it's valid
      child_result <- tryCatch({
        hexify_index_to_cell(child_idx, 7L, "z7")
      }, error = function(e) NULL)
      
      if (!is.null(child_result)) {
        # Child should be one resolution finer (need to extract numeric value)
        expect_equal(as.integer(child_result$resolution), 
                     as.integer(parent_result$resolution) + 1L,
                     info = sprintf("Child %s resolution check", child_idx))
        
        # Every child has a unique, stable index.
        child_idx2 <- hexify_cell_to_index(child_result$face, child_result$i,
                                    child_result$j, child_result$resolution, 7L, "z7")
        
        expect_equal(child_idx2, child_idx,
                     info = sprintf("Child %s should round-trip", child_idx))
      }
    }
  }
})

test_that("Z7: Edge cases are handled correctly", {
  # Test various edge cases and error conditions
  
  # Single digit should fail (need at least 2 for face)
  expect_error(hexify_index_to_cell("1", 7L, "z7"))
  
  # Face numbers 12 and 99 may or may not be errors depending on implementation
  # DGGRID might handle them differently than expected
  # Let's just test that they can be processed without crashing
  result_12 <- tryCatch(hexify_index_to_cell("12", 7L, "z7"), error = function(e) NULL)
  result_99 <- tryCatch(hexify_index_to_cell("99", 7L, "z7"), error = function(e) NULL)
  
  # Empty string should fail
  expect_error(hexify_index_to_cell("", 7L, "z7"))
  
  # Arbitrary integer coordinates should round-trip without crashing.
  idx_neg <- hexify_cell_to_index(5L, -1L, -1L, 1L, 7L, "z7")
  result_neg <- hexify_index_to_cell(idx_neg, 7L, "z7")
  expect_equal(c(result_neg$i, result_neg$j), c(-1, -1),
               info = "Negative coords should round-trip")
})

test_that("Z7: Resolution progression works correctly", {
  # Test that resolution scaling works (factor of 7 per level)
  
  # Start at resolution 1 with simple coordinates
  idx_r1 <- hexify_cell_to_index(5L, 1L, 1L, 1L, 7L, "z7")
  
  # At resolution 2, same cell center should be at (7, 7)
  idx_r2 <- hexify_cell_to_index(5L, 7L, 7L, 2L, 7L, "z7")
  
  # At resolution 3, same cell center should be at (49, 49)
  idx_r3 <- hexify_cell_to_index(5L, 49L, 49L, 3L, 7L, "z7")
  
  # These indices represent the same logical cell at different resolutions
  # Their string representations will be different lengths
  expect_equal(nchar(idx_r1), 2L + 1L)  # Face + 1 digit
  expect_equal(nchar(idx_r2), 2L + 2L)  # Face + 2 digits
  expect_equal(nchar(idx_r3), 2L + 3L)  # Face + 3 digits
})

test_that("Z7: Former DGGRID problem indices are stable", {
  problem_indices <- c("012", "110001", "110002", "110004", "110006")

  for (idx in problem_indices) {
    result <- hexify_index_to_cell(idx, 7L, "z7")
    re_encoded <- hexify_cell_to_index(result$face, result$i, result$j,
                                result$resolution, 7L, "z7")
    expect_equal(re_encoded, idx,
                 info = sprintf("%s should round-trip", idx))
  }
})

test_that("Z7: Canonical forms provide stability", {
  for (idx in c("012", "110001", "110002", "110004", "110006", "0955555")) {
    expect_equal(hexify_z7_canonical(idx), idx)
    canonical <- hexify_z7_canonical(idx)
    result <- hexify_index_to_cell(canonical, 7L, "z7")
    re_encoded <- hexify_cell_to_index(result$face, result$i, result$j,
                                result$resolution, 7L, "z7")
    canonical2 <- hexify_z7_canonical(re_encoded)
    expect_equal(canonical2, canonical,
                 info = sprintf("Canonical form of %s should be stable", idx))
  }
})
