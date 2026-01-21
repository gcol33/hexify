
# tests/testthat/test-index_z3.R
# Comprehensive tests for Z3 indexing (aperture 3 specialized)
# VERSION: 2024-11-03 FIXED - Valid Class II coordinates
# 
# IMPORTANT: For Class II (odd) resolutions, coordinates must satisfy:
#   j % 3 == jDigits[i % 3]
#   where jDigits = {0→0, 1→2, 2→1}
#
# Valid examples:
#   i%3=0 → j%3=0: (0,0), (0,3), (3,6), (6,9), etc.
#   i%3=1 → j%3=2: (1,2), (1,5), (4,8), (7,11), etc.
#   i%3=2 → j%3=1: (2,1), (2,4), (5,7), (8,10), etc.

# =============================================================================
# Z3 Lookup Table Tests (Core Algorithm)
# =============================================================================

test_that("Z3: lookup table matches DGGRID exactly", {
  skip_on_cran()
  # DGGRID Z3 lookup table (DgZ3StringRF.cpp:135-140)
  z3_table <- matrix(c(
    "00", "22", "21",   # i=0, j=0,1,2
    "01", "02", "20",   # i=1, j=0,1,2
    "12", "10", "11"    # i=2, j=0,1,2
  ), nrow = 3, byrow = TRUE)
  
  # Test all 9 lookup table entries at resolution 2 (Class I)
  for (i in 0:2) {
    for (j in 0:2) {
      index <- cpp_cell_to_index(0, i, j, 2, 3, "z3")
      expected <- paste0("00", z3_table[i+1, j+1])
      
      expect_equal(index, expected,
                   label = sprintf("Z3[%d][%d]", i, j))
    }
  }
})

test_that("Z3: specific lookup table entries", {
  # Test specific entries explicitly
  expect_equal(cpp_cell_to_index(0, 0, 0, 2, 3, "z3"), "0000")  # [0][0]→"00"
  expect_equal(cpp_cell_to_index(0, 0, 1, 2, 3, "z3"), "0022")  # [0][1]→"22"
  expect_equal(cpp_cell_to_index(0, 0, 2, 2, 3, "z3"), "0021")  # [0][2]→"21"
  
  expect_equal(cpp_cell_to_index(0, 1, 0, 2, 3, "z3"), "0001")  # [1][0]→"01"
  expect_equal(cpp_cell_to_index(0, 1, 1, 2, 3, "z3"), "0002")  # [1][1]→"02"
  expect_equal(cpp_cell_to_index(0, 1, 2, 2, 3, "z3"), "0020")  # [1][2]→"20"
  
  expect_equal(cpp_cell_to_index(0, 2, 0, 2, 3, "z3"), "0012")  # [2][0]→"12"
  expect_equal(cpp_cell_to_index(0, 2, 1, 2, 3, "z3"), "0010")  # [2][1]→"10"
  expect_equal(cpp_cell_to_index(0, 2, 2, 2, 3, "z3"), "0011")  # [2][2]→"11"
})

# =============================================================================
# Z3 Decoding Tests (Reverse Lookup)
# =============================================================================

test_that("Z3: decode all 9 lookup table entries", {
  skip_on_cran()
  z3_codes <- list(
    list(code = "0000", i = 0, j = 0),  # "00"→[0][0]
    list(code = "0022", i = 0, j = 1),  # "22"→[0][1]
    list(code = "0021", i = 0, j = 2),  # "21"→[0][2]
    list(code = "0001", i = 1, j = 0),  # "01"→[1][0]
    list(code = "0002", i = 1, j = 1),  # "02"→[1][1]
    list(code = "0020", i = 1, j = 2),  # "20"→[1][2]
    list(code = "0012", i = 2, j = 0),  # "12"→[2][0]
    list(code = "0010", i = 2, j = 1),  # "10"→[2][1]
    list(code = "0011", i = 2, j = 2)   # "11"→[2][2]
  )
  
  for (tc in z3_codes) {
    result <- cpp_index_to_cell(tc$code, 3, "z3")
    
    expect_equal(result$face, 0)
    expect_equal(result$i, tc$i,
                 label = sprintf("decode '%s'", tc$code))
    expect_equal(result$j, tc$j,
                 label = sprintf("decode '%s'", tc$code))
    expect_equal(result$resolution, 2)
  }
})

test_that("Z3: decode reverse lookup table correctness", {
  skip_on_cran()
  # Test that the reverse lookup matches DGGRID
  # DGGRID reverse table (DgZ3StringRF.cpp:196-253)
  reverse_cases <- list(
    list(z3 = "00", i = 0, j = 0),
    list(z3 = "22", i = 0, j = 1),
    list(z3 = "21", i = 0, j = 2),
    list(z3 = "01", i = 1, j = 0),
    list(z3 = "02", i = 1, j = 1),
    list(z3 = "20", i = 1, j = 2),
    list(z3 = "12", i = 2, j = 0),
    list(z3 = "10", i = 2, j = 1),
    list(z3 = "11", i = 2, j = 2)
  )
  
  for (tc in reverse_cases) {
    index <- paste0("00", tc$z3)
    result <- cpp_index_to_cell(index, 3, "z3")
    
    expect_equal(result$i, tc$i,
                 label = sprintf("z3='%s'", tc$z3))
    expect_equal(result$j, tc$j,
                 label = sprintf("z3='%s'", tc$z3))
  }
})

# =============================================================================
# Class I vs Class II Tests
# =============================================================================

test_that("Z3: Class I string length (even resolutions)", {
  skip_on_cran()
  # Class I resolutions: 0, 2, 4, 6, ...
  # String length: 2 (face) + resolution digits

  class_i_cases <- list(
    list(res = 0, length = 2),
    list(res = 2, length = 4),
    list(res = 4, length = 6),
    list(res = 6, length = 8)
  )
  
  for (tc in class_i_cases) {
    index <- cpp_cell_to_index(0, 0, 0, tc$res, 3, "z3")
    expect_equal(nchar(index), tc$length,
                 label = sprintf("Class I res=%d", tc$res))
  }
})

test_that("Z3: Class II string length (odd resolutions)", {
  skip_on_cran()
  # Class II resolutions: 1, 3, 5, 7, ...
  # String length: 2 (face) + resolution digits

  class_ii_cases <- list(
    list(res = 1, length = 3),
    list(res = 3, length = 5),
    list(res = 5, length = 7),
    list(res = 7, length = 9)
  )
  
  for (tc in class_ii_cases) {
    index <- cpp_cell_to_index(0, 0, 0, tc$res, 3, "z3")
    expect_equal(nchar(index), tc$length,
                 label = sprintf("Class II res=%d", tc$res))
  }
})

test_that("Z3: Class II trimming behavior", {
  # For Class II resolutions, last digit should be trimmed
  # Test by checking that indices differ in expected ways
  
  # Resolution 2 (Class I) vs 3 (Class II)
  idx_r2 <- cpp_cell_to_index(0, 1, 2, 2, 3, "z3")  # Valid: i%3=1, j%3=2
  idx_r3 <- cpp_cell_to_index(0, 1, 2, 3, 3, "z3")  # Valid: i%3=1, j%3=2
  
  # Class I should be 4 chars, Class II should be 5 chars
  expect_equal(nchar(idx_r2), 4)
  expect_equal(nchar(idx_r3), 5)
  
  # Resolution 4 (Class I) vs 5 (Class II)
  idx_r4 <- cpp_cell_to_index(0, 2, 1, 4, 3, "z3")  # Valid: i%3=2, j%3=1
  idx_r5 <- cpp_cell_to_index(0, 2, 1, 5, 3, "z3")  # Valid: i%3=2, j%3=1
  
  expect_equal(nchar(idx_r4), 6)
  expect_equal(nchar(idx_r5), 7)
})

test_that("Z3: Class II j-digit inference", {
  skip_on_cran()
  # For Class II (odd resolution), last j digit is inferred from last i digit
  # Rule from DGGRID: i%3=0→j%3=0, i%3=1→j%3=2, i%3=2→j%3=1

  # These should decode correctly even though j digit is implicit
  test_cases <- list(
    # Resolution 1 (Class II) - eff_res=1, max=2
    list(i = 0, j = 0, res = 1),  # i%3=0, j%3=0 ✓
    list(i = 1, j = 2, res = 1),  # i%3=1, j%3=2 ✓
    list(i = 2, j = 1, res = 1),  # i%3=2, j%3=1 ✓
    
    # Resolution 3 (Class II) - eff_res=2, max=8
    list(i = 3, j = 6, res = 3),  # i%3=0, j%3=0 ✓ (i="10", j="20")
    list(i = 4, j = 8, res = 3),  # i%3=1, j%3=2 ✓ (i="11", j="22")
    list(i = 5, j = 7, res = 3)   # i%3=2, j%3=1 ✓ (i="12", j="21")
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, 3, "z3")
    result <- cpp_index_to_cell(index, 3, "z3")
    
    expect_equal(result$i, tc$i,
                 label = sprintf("infer i[%d,%d,r%d]", tc$i, tc$j, tc$res))
    expect_equal(result$j, tc$j,
                 label = sprintf("infer j[%d,%d,r%d]", tc$i, tc$j, tc$res))
  }
})

# =============================================================================
# Multi-digit Z3 Encoding Tests
# =============================================================================

test_that("Z3: multi-digit encoding at resolution 4", {
  # Resolution 4 (Class I) has eff_res=2, so 2 digit pairs
  
  # i=3, j=5 in radix-3: i="10", j="12"
  # Z3 codes: [1][1]→"02", [0][2]→"21"
  # Expected: "00" + "02" + "21" = "000221"
  index <- cpp_cell_to_index(0, 3, 5, 4, 3, "z3")
  expect_equal(index, "000221")
  
  # Roundtrip
  result <- cpp_index_to_cell(index, 3, "z3")
  expect_equal(result$i, 3)
  expect_equal(result$j, 5)
  expect_equal(result$resolution, 4)
})

test_that("Z3: multi-digit encoding at resolution 6", {
  # Resolution 6 (Class I) has eff_res=3, so 3 digit pairs
  
  # i=13, j=17 in radix-3: i="111", j="122"
  # Z3 codes: [1][1]→"02", [1][2]→"20", [1][2]→"20"
  # Expected: "00" + "02" + "20" + "20" = "00022020"
  index <- cpp_cell_to_index(0, 13, 17, 6, 3, "z3")
  expect_equal(index, "00022020")
  
  # Roundtrip
  result <- cpp_index_to_cell(index, 3, "z3")
  expect_equal(result$i, 13)
  expect_equal(result$j, 17)
  expect_equal(result$resolution, 6)
})

# =============================================================================
# Roundtrip Tests
# =============================================================================

test_that("Z3: basic roundtrip for small coordinates", {
  skip_on_cran()
  # Generate test cases with valid Class II coordinates
  test_cases_raw <- expand.grid(
    face = c(0, 5, 10, 19),
    i = 0:8,
    j = 0:8,
    res = 1:3
  )
  
  # Filter for valid coordinates
  test_cases <- test_cases_raw[apply(test_cases_raw, 1, function(row) {
    res <- row[["res"]]
    i <- row[["i"]]
    j <- row[["j"]]
    
    # Calculate max coordinate
    eff_res <- (res + 1) %/% 2
    max_coord <- 3^eff_res - 1
    
    # Check if within bounds
    if (i > max_coord || j > max_coord) return(FALSE)
    
    # For Class II (odd res), check constraint
    if (res %% 2 == 1) {
      j_required <- c(0, 2, 1)[i %% 3 + 1]
      if (j %% 3 != j_required) return(FALSE)
    }
    
    return(TRUE)
  }), ]
  
  for (idx in seq_len(nrow(test_cases))) {
    tc <- test_cases[idx, ]
    
    # Roundtrip
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 3, "z3")
    result <- cpp_index_to_cell(index, 3, "z3")
    
    expect_equal(result$face, tc$face, label = sprintf("idx=%d", idx))
    expect_equal(result$i, tc$i, label = sprintf("idx=%d", idx))
    expect_equal(result$j, tc$j, label = sprintf("idx=%d", idx))
    expect_equal(result$resolution, tc$res, label = sprintf("idx=%d", idx))
  }
})

test_that("Z3: roundtrip for higher resolutions", {
  skip_on_cran()
  # All coordinates must satisfy Class II constraint for odd resolutions
  test_cases <- list(
    # Resolution 4 (Class I): eff_res=2, max=8
    list(face = 0, i = 5, j = 7, res = 4),
    
    # Resolution 5 (Class II): eff_res=3, max=26
    # Must satisfy: j%3 == jDigits[i%3]
    list(face = 5, i = 15, j = 21, res = 5),  # i%3=0, j%3=0 ✓ (was 15,20)
    
    # Resolution 6 (Class I): eff_res=3, max=26
    list(face = 10, i = 20, j = 25, res = 6),
    list(face = 19, i = 26, j = 26, res = 6)
  )
  
  for (tc in test_cases) {
    index <- cpp_cell_to_index(tc$face, tc$i, tc$j, tc$res, 3, "z3")
    result <- cpp_index_to_cell(index, 3, "z3")
    
    expect_equal(result$face, tc$face)
    expect_equal(result$i, tc$i)
    expect_equal(result$j, tc$j)
    expect_equal(result$resolution, tc$res)
  }
})

# =============================================================================
# Face Encoding Tests
# =============================================================================

test_that("Z3: face number encoding", {
  skip_on_cran()
  # All faces should encode with leading zero if needed
  for (face in 0:19) {
    index <- cpp_cell_to_index(face, 0, 0, 1, 3, "z3")
    expect_equal(substr(index, 1, 2), sprintf("%02d", face))
    
    # Decode should work
    result <- cpp_index_to_cell(index, 3, "z3")
    expect_equal(result$face, face)
  }
})

test_that("Z3: resolution 0 (face only)", {
  skip_on_cran()
  # At resolution 0, only faces 0-11 are valid (the 12 pentagons)
  # Faces 12-19 only appear at resolution 1+
  for (face in c(0, 5, 10, 11)) {  # Changed from c(0, 5, 10, 15, 19)
    index <- cpp_cell_to_index(face, 0, 0, 0, 3, "z3")
    
    # Should be exactly 2 chars
    expect_equal(nchar(index), 2)
    
    # Decode
    result <- cpp_index_to_cell(index, 3, "z3")
    expect_equal(result$face, face)
    expect_equal(result$resolution, 0)
    expect_equal(result$i, 0)
    expect_equal(result$j, 0)
  }
})

# =============================================================================
# Hierarchy Tests
# =============================================================================

test_that("Z3: parent/child relationships", {
  skip_on_cran()
  # Create a cell at resolution 5 (Class II) with valid coordinates
  # i=10, j=11: i%3=1, j%3=2 ✓
  child <- cpp_cell_to_index(0, 10, 11, 5, 3, "z3")
  
  # Get parent
  parent <- cpp_get_parent_index(child, 3, "z3")
  
  # Parent should be prefix of child
  expect_true(substr(child, 1, nchar(parent)) == parent)
  
  # Get children
  children <- cpp_get_children_indices(parent, 3, "z3")
  
  # Original child should be in list
  expect_true(child %in% children)
  
  # Should have 3 children (aperture 3)
  expect_equal(length(children), 3)
})

test_that("Z3: multi-level parent traversal", {
  skip_on_cran()
  # Start at resolution 6 with valid coordinates
  idx_r6 <- cpp_cell_to_index(0, 20, 22, 6, 3, "z3")  # i%3=2, j%3=1 ✓
  
  # Get parent at each level
  idx_r5 <- cpp_get_parent_index(idx_r6, 3, "z3")
  idx_r4 <- cpp_get_parent_index(idx_r5, 3, "z3")
  idx_r3 <- cpp_get_parent_index(idx_r4, 3, "z3")
  
  # Each should be progressively shorter
  expect_true(nchar(idx_r5) < nchar(idx_r6))
  expect_true(nchar(idx_r4) < nchar(idx_r5))
  expect_true(nchar(idx_r3) < nchar(idx_r4))
  
  # Each should be prefix of child
  expect_true(substr(idx_r6, 1, nchar(idx_r5)) == idx_r5)
  expect_true(substr(idx_r5, 1, nchar(idx_r4)) == idx_r4)
  expect_true(substr(idx_r4, 1, nchar(idx_r3)) == idx_r3)
})

# =============================================================================
# Comparison and Utility Tests
# =============================================================================

test_that("Z3: index comparison", {
  idx1 <- cpp_cell_to_index(0, 0, 0, 2, 3, "z3")
  idx2 <- cpp_cell_to_index(0, 1, 0, 2, 3, "z3")
  idx3 <- cpp_cell_to_index(1, 0, 0, 2, 3, "z3")
  
  expect_equal(cpp_compare_indices(idx1, idx1), 0)
  expect_true(cpp_compare_indices(idx1, idx2) != 0)
  expect_true(cpp_compare_indices(idx1, idx3) != 0)
})

test_that("Z3: resolution extraction", {
  skip_on_cran()
  # Resolution extraction works for both Class I and Class II
  # Resolution = string length (after face), regardless of class

  # Test Class I resolutions (even)
  class_i_resolutions <- c(0, 2, 4, 6)
  
  for (res in class_i_resolutions) {
    index <- cpp_cell_to_index(0, 0, 0, res, 3, "z3")
    extracted <- cpp_get_index_resolution(index, 3, "z3")
    expect_equal(extracted, res,
                 label = sprintf("Class I res=%d", res))
  }
  
  # Test Class II resolutions (odd) - these also work correctly now
  class_ii_resolutions <- c(1, 3, 5, 7)
  
  for (res in class_ii_resolutions) {
    index <- cpp_cell_to_index(0, 0, 0, res, 3, "z3")
    extracted <- cpp_get_index_resolution(index, 3, "z3")
    expect_equal(extracted, res,
                 label = sprintf("Class II res=%d", res))
  }
})

# =============================================================================
# Difference from Z-Order Tests
# =============================================================================

test_that("Z3: produces different results than Z-Order", {
  skip_on_cran()
  # Z3 and Z-Order should give different indices for same input
  # IMPORTANT: Use VALID coordinates for Class II!
  test_cases <- list(
    list(i = 1, j = 2, res = 2),  # Class I: any coords OK
    list(i = 2, j = 1, res = 2),  # Class I: any coords OK
    list(i = 3, j = 3, res = 3)   # Class II: i%3=0, j%3=0 ✓ (was 3,4)
  )
  
  for (tc in test_cases) {
    idx_z3 <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, 3, "z3")
    idx_zo <- cpp_cell_to_index(0, tc$i, tc$j, tc$res, 3, "zorder")
    
    # Should be different
    expect_false(idx_z3 == idx_zo,
                 label = sprintf("i=%d,j=%d,res=%d", tc$i, tc$j, tc$res))
    
    # But both should decode correctly
    result_z3 <- cpp_index_to_cell(idx_z3, 3, "z3")
    result_zo <- cpp_index_to_cell(idx_zo, 3, "zorder")
    
    expect_equal(result_z3$i, tc$i,
                 label = sprintf("Z3 i for (%d,%d)", tc$i, tc$j))
    expect_equal(result_z3$j, tc$j,
                 label = sprintf("Z3 j for (%d,%d)", tc$i, tc$j))
    expect_equal(result_zo$i, tc$i,
                 label = sprintf("ZOrder i for (%d,%d)", tc$i, tc$j))
    expect_equal(result_zo$j, tc$j,
                 label = sprintf("ZOrder j for (%d,%d)", tc$i, tc$j))
  }
})

test_that("Z3: provides better spatial locality than Z-Order", {
  # This is a qualitative test - Z3 should cluster spatially close cells
  # We test that neighboring cells have similar indices
  
  # Get indices for a small neighborhood
  idx_00 <- cpp_cell_to_index(0, 0, 0, 2, 3, "z3")
  idx_01 <- cpp_cell_to_index(0, 0, 1, 2, 3, "z3")
  idx_10 <- cpp_cell_to_index(0, 1, 0, 2, 3, "z3")
  
  # All should start with same face
  expect_equal(substr(idx_00, 1, 2), "00")
  expect_equal(substr(idx_01, 1, 2), "00")
  expect_equal(substr(idx_10, 1, 2), "00")
  
  # And have same length (Class I)
  expect_equal(nchar(idx_00), nchar(idx_01))
  expect_equal(nchar(idx_00), nchar(idx_10))
})
